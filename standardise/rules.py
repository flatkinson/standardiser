"""Apply standardisation rules"""

########################################################################

# Module imports...

from __future__ import print_function

import logging; logger = logging.getLogger(__name__)

import os
import re
from copy import deepcopy

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

rules_fh = open(os.path.join(os.path.dirname(__file__), data_dir_name, rules_file_name))

rule_set = [{"n": n, "description": x[1], "SMARTS": x[0], "rxn": AllChem.ReactionFromSmarts(x[0])} for n, x in enumerate([x[0] for x in [re.findall("^\s*([^#]\S+)\t+(.*?)\s*$", line) for line in rules_fh] if x], 1)]

########################################################################

def setAllHsExplicit(mol):

    for atom in mol.GetAtoms():

        atom.SetNumExplicitHs(atom.GetTotalNumHs())

        atom.SetNoImplicit(True)

# setAllHsExplicit

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

			logger.debug("SanitizeMol failed for rule {} ({}): '{}'".format(rule['n'], rule['description'], e.message.rstrip()))

			break
	
		mol = new_mol

                n_hits_for_rule += 1

            if n_hits_for_rule:
                
                logger.debug("rule {} ({}) applied {} time(s) on pass {}".format(rule['n'], rule['description'], n_hits_for_rule, n_pass))

                n_hits_for_pass += n_hits_for_rule

        logger.debug("...total of {} applications in pass: {}".format(n_hits_for_pass, "will continue..." if n_hits_for_pass else "finished."))

        if n_hits_for_pass == 0: break

    setAllHsExplicit(mol)

    return mol

# apply_old

######

def apply_rule(mol, rule, verbose=False):

    if verbose: logger.debug("apply_rule> applying rule {} ({})...".format(rule["n"], rule["description"]))

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
        
    logging.debug("apply_rule {} ({})> maximum number of passes reached; current number of mols is {}".format(rule["n"], rule["description"], len(mols)))
    
    return mols[0]

# apply_rule

######

def apply(mol):

    for n_pass in range(1, max_passes+1):

        logger.debug("apply> starting pass {}...".format(n_pass))

        n_hits_for_pass = 0

        for rule in rule_set:

	    product = apply_rule(mol, rule)

	    if product:

                logger.debug("apply> hit for rule {} ({})...".format(rule["n"], rule["description"]))

                mol = product

                n_hits_for_pass += 1

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

            logger.info("rule {} ('{}') applied".format(rule['n'], rule['description']))

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

    coord_map = {new_idx: rdGeometry.Point2D(conf.GetAtomPosition(old_idx).x, conf.GetAtomPosition(old_idx).y) for old_idx, new_idx in zip(old_match, new_match)}

    AllChem.Compute2DCoords(new_mol, clearConfs=True, coordMap=coord_map, canonOrient=False)

    return old_mol, old_match, new_mol, new_match

# demo

########################################################################
# End
########################################################################
