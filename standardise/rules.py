"""Apply standardisation rules"""

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

    for atom in mol.GetAtoms():

        atom.SetNumExplicitHs(atom.GetTotalNumHs())

        atom.SetNoImplicit(True)

# setAllHsExplicit

######

def apply_rule(mol, rule, verbose=False):

    if verbose: logger.debug("apply_rule> applying rule {} '{}'...".format(rule["n"], rule["name"]))

    mols = [mol]

    changed = False

    for n_pass in range(1, max_passes+1):

        if verbose: logging.debug("apply_rule> starting pass {}...".format(n_pass))

        products = {}

        for mol in mols:

            for product in [x[0] for x in rule["rxn"].RunReactants((mol,))]:

                try:
                    Chem.SanitizeMol(product)
                except ValueError as error:
                    if error.message.startswith("Sanitization error"): continue

                products.setdefault(Chem.MolToSmiles(product, isomericSmiles=True), product)

        if products:

            changed = True

            if (verbose): logging.debug("apply_rule> there are {} products: will continue".format(len(products.values())))

            mols = products.values()

        else:

            if (verbose): logging.debug("apply_rule> there were no products: will return")

            return mols[0] if changed else None

    logging.debug("apply_rule {} '{}'> maximum number of passes reached; current number of mols is {}".format(rule["n"], rule["name"], len(mols)))

    return mols[0]

# apply_rule

######

def apply(mol, first_only=False, verbose=False, output_rules_applied=None):

    logger.debug("apply> mol = '{}'".format(Chem.MolToSmiles(mol)))

    rules_applied = []

    for n_pass in range(1, max_passes+1):

        logger.debug("apply> starting pass {}...".format(n_pass))

        n_hits_for_pass = 0

        for rule in rule_set:

            product = apply_rule(mol, rule, verbose)

            if product:

                logger.info("rule {} '{}' applied on pass {}".format(rule["n"], rule["name"], n_pass))

                mol = product

                if output_rules_applied is not None: output_rules_applied.append(rule["n"])

                if first_only: break

                n_hits_for_pass += 1

        if product and first_only: break

        logger.debug("...total of {} hits in pass: {}".format(n_hits_for_pass, "will continue..." if n_hits_for_pass else "finished."))

        if n_hits_for_pass == 0: break

    setAllHsExplicit(mol)

    return mol

# apply

######

def demo(old_mol):

    new_mol = None

    for rule in rule_set:

        products = rule["rxn"].RunReactants((old_mol,))

        if len(products):

            new_mol = products[0][0]

            logger.info("rule {} '{}' applied".format(rule['n'], rule['name']))

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

######

def apply_old(mol):

    for n_pass in range(1, max_passes+1):

        logger.debug("starting pass {}...".format(n_pass))

        n_hits_for_pass = 0

        for rule in rule_set:

            n_hits_for_rule = 0

            while True:

                products = rule["rxn"].RunReactants((mol,))

                if not len(products): break

                new_mol = products[0][0]

                try:
                    Chem.SanitizeMol(new_mol)

                except ValueError as e:

                    logger.debug("SanitizeMol failed for rule {} '{}': '{}'".format(rule['n'], rule['name'], e.message.rstrip()))

                    break

                mol = new_mol

                n_hits_for_rule += 1

            if n_hits_for_rule:

                logger.debug("rule {} '{}' applied {} time(s) on pass {}".format(rule['n'], rule['name'], n_hits_for_rule, n_pass))

                n_hits_for_pass += n_hits_for_rule

        logger.debug("...total of {} applications in pass: {}".format(n_hits_for_pass, "will continue..." if n_hits_for_pass else "finished."))

        if n_hits_for_pass == 0: break

    setAllHsExplicit(mol)

    return mol

# apply_old

########################################################################
# End
########################################################################
