####################################################################################################
# 
# Copyright [2014] EMBL - European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in
# compliance with the License.  You may obtain a copy of
# the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software distributed under the
# License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
# either express or implied. See the License for the specific language governing permissions
# and limitations under the License.
# 
####################################################################################################

"""
Module to apply rule-based standardisations.
"""

########################################################################

# Module imports...

from __future__ import print_function

import logging
logger = logging.getLogger(__name__)

import os
import re
import csv
from itertools import ifilterfalse

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import rdGeometry

from standardiser.utils import StandardiseException

########################################################################

# Module configuration...

rules_file_name = "rules.dat"

data_dir_name = "data"

max_passes = 10

########################################################################

# Module initialization...

with open(os.path.join(os.path.dirname(__file__), data_dir_name, rules_file_name)) as rules_file:

    reader = csv.reader(ifilterfalse(lambda x: re.match(r"^\s*(?:#|$)", x), rules_file), delimiter="\t") # SMARTS and name, tab-seperated

    rule_set = [{"n": n, "SMARTS": x[0], "rxn": AllChem.ReactionFromSmarts(x[0]), "name": x[1]} for n, x in enumerate(reader, 1)]

########################################################################

def setAllHsExplicit(mol):

    """
    Make all hydrogens explicit.
   
    #TODO: In the earliest versions of this project, this hack was necessary to get things working as I expected.
    That was several versions of RDKit ago, however, and I really need to check whether this is still necessary.
    """

    for atom in mol.GetAtoms():

        atom.SetNumExplicitHs(atom.GetTotalNumHs())

        atom.SetNoImplicit(True)

# setAllHsExplicit

######

def apply_rule(mol, rule, verbose=False):

    """
    Apply a single rule to the input molecule.

    Please see the IPython Notebook 'issue_01' for an explanation of why things are done this way.
    """

    if verbose: logger.debug("apply_rule> applying rule {n} '{name}'...".format(n=rule["n"], name=rule["name"]))

    mols = [mol]

    changed = False

    for n_pass in range(1, max_passes+1):

        if verbose: logging.debug("apply_rule> starting pass {n}...".format(n=n_pass))

        products = {}

        for mol in mols:

            for product in [x[0] for x in rule["rxn"].RunReactants((mol,))]:

                try:

                    Chem.SanitizeMol(product)

                    smiles = Chem.MolToSmiles(product, isomericSmiles=True)

                except ValueError as error:

                    continue # We are assuming this simply means an unphysical molecule has been generated

                if smiles in products: continue # Keep only new structures
                
                products[smiles] = product

        if products:

            changed = True

            if (verbose): logging.debug("apply_rule> there are {n} products: will continue".format(n=len(products.values())))

            mols = products.values() # Update list of mols

        else:

            if (verbose): logging.debug("apply_rule> there were no products: will return")

            return mols[0] if changed else None

    logging.debug("apply_rule {n} '{name}'> maximum number of passes reached; current number of mols is {m}".format(n=rule["n"], name=rule["name"], m=len(mols)))

    return mols[0]

# apply_rule

######

def apply(mol, first_only=False, verbose=False, output_rules_applied=None):

    """
    Apply all rules to the input molecule.
    """

    logger.debug("apply> mol = '{smi}'".format(smi=Chem.MolToSmiles(mol)))

    rules_applied = []

    for n_pass in range(1, max_passes+1):

        logger.debug("apply> starting pass {n}...".format(n=n_pass))

        n_hits_for_pass = 0

        for rule in rule_set:

            product = apply_rule(mol, rule, verbose)

            if product:

                logger.debug("rule {n} '{name}' applied on pass {m}".format(n=rule["n"], name=rule["name"], m=n_pass))

                mol = product

                if output_rules_applied is not None: output_rules_applied.append(rule["n"])

                if first_only: break

                n_hits_for_pass += 1

        if product and first_only: break

        logger.debug("...total of {n} hits in pass: {m}".format(n=n_hits_for_pass, m="will continue..." if n_hits_for_pass else "finished."))

        if n_hits_for_pass == 0: break

    setAllHsExplicit(mol) 

    return mol

# apply

######

def demo(old_mol):

    """
    Utility function for illustrating the application of rules.
    """

    new_mol = None

    for rule in rule_set:

        products = rule["rxn"].RunReactants((old_mol,))

        if len(products):

            new_mol = products[0][0]

            logger.debug("rule {n} '{name}' applied".format(n=rule['n'], name=rule['name']))

            break

    if not new_mol:

        logger.warn("No hits for mol!")

        return None, None, None, None

    Chem.SanitizeMol(new_mol)

    AllChem.Compute2DCoords(old_mol)

    conf = old_mol.GetConformer()

    old_pat, new_pat = rule["SMARTS"].split(">>")

    old_match = old_mol.GetSubstructMatch(Chem.MolFromSmarts(old_pat))
    new_match = new_mol.GetSubstructMatch(Chem.MolFromSmarts(new_pat))

    #@ coord_map = {new_idx: rdGeometry.Point2D(conf.GetAtomPosition(old_idx).x, conf.GetAtomPosition(old_idx).y) for old_idx, new_idx in zip(old_match, new_match)}
    coord_map = dict((new_idx, rdGeometry.Point2D(conf.GetAtomPosition(old_idx).x, conf.GetAtomPosition(old_idx).y)) for old_idx, new_idx in zip(old_match, new_match))  # python2.6-compatible

    AllChem.Compute2DCoords(new_mol, clearConfs=True, coordMap=coord_map, canonOrient=False)

    return old_mol, old_match, new_mol, new_match

# demo

########################################################################
# End
########################################################################
