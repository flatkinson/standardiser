"""Apply standardisation procedure"""

########################################################################

from __future__ import print_function

import logging
logger = logging.getLogger(__name__)

import rdkit
from rdkit import Chem

from . import break_bonds, neutralise, rules, unsalt

########################################################################
#
# Module configuration...
#

errors = {
    "not_built":       "RDKit could not build mol",
    "no_non_salt":     "No non-salt/solvate components",
    "multi_component": "Multiple non-salt/solvate components"
}

########################################################################
#
# Module initialization...
#

########################################################################

class StandardiseException(Exception):

    def __init__(self, name):

        self.name = name

        self.message = errors[name]

        self.args = (errors[name], )

# StandardiseException

######

def apply(input_mol, output_rules_applied=None): #@@@

    # Get input molecule...

    if type(input_mol) == rdkit.Chem.rdchem.Mol:

        mol = input_mol

        input_type = 'mol'

    else:

        mol = Chem.MolFromMolBlock(input_mol)

        if not mol:

            mol = Chem.MolFromSmiles(input_mol)

            if not mol:

                raise StandardiseException("not_built")

            else:

                input_type = 'smi'
        else:

            input_type = 'sdf'

    ######

    # Get disconnected fragments...

    non_salt_frags = []

    mol = break_bonds.apply(mol)

    for n, frag in enumerate(Chem.GetMolFrags(mol, asMols=True), 1):

        logger.debug("Starting fragment {n} '{smi}'...".format(n=n, smi=Chem.MolToSmiles(frag)))

        logger.debug("1) Check for non-organic elements...")

        if unsalt.is_nonorganic(frag): continue

        logger.debug("2) Attempting to neutralise (first pass)...")

        frag = neutralise.apply(frag)

        logger.debug("3) Applying rules...")

        frag = rules.apply(frag, output_rules_applied=output_rules_applied) #@@@

        logger.debug("4) Attempting to neutralise (second pass)...")

        frag = neutralise.apply(frag)

        logger.debug("5) Checking if frag is a salt/solvate...")

        if unsalt.is_salt(frag): continue

        logger.debug("...fragment kept.")

        non_salt_frags.append(frag)

    if len(non_salt_frags) == 0:

        raise StandardiseException("no_non_salt")

    if len(non_salt_frags) > 1:

        raise StandardiseException("multi_component")

    parent = non_salt_frags[0]

    ######

    # Return parent in same format as input...

    if input_type == 'mol':

        return parent

    elif input_type == 'sdf':

        return Chem.MolToMolBlock(parent)

    else:  # input_type == 'smi'

        return Chem.MolToSmiles(parent, isomericSmiles=True)

# apply

########################################################################
# End
########################################################################
