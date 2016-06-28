################################################################################################################################
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
################################################################################################################################

"""
Neutralise charges
"""

####################################################################################################

import make_logger
logger = make_logger.run(__name__)

import copy

from rdkit import Chem

from .utils import StandardiseException, sanity_check

####################################################################################################

# Module configuration...

pos_smarts    = "[+!H0!$(*~[-])]"
quat_smarts   = "[+H0!$(*~[-])]"
neg_smarts    = "[-!$(*~[+H0])]"
acid_smarts   = "[$([O-][C,P,S]=O),$([n-]1nnnc1),$(n1[n-]nnc1)]"
acid_h_smarts = "[$([OH][C,P,S]=O),$([n-]1nnnc1),$(n1[n-]nnc1)]"

####################################################################################################

# Module initialization...

pos_pat    = Chem.MolFromSmarts(pos_smarts)
quat_pat   = Chem.MolFromSmarts(quat_smarts)
neg_pat    = Chem.MolFromSmarts(neg_smarts)
acid_pat   = Chem.MolFromSmarts(acid_smarts)
acid_h_pat = Chem.MolFromSmarts(acid_h_smarts)

####################################################################################################

def formal_charge(mol):

    return sum(x.GetFormalCharge() for x in mol.GetAtoms())

# def formal_charge

######

# NB This solved a problem in early iterations: newer RDKit releases may have made it unnecessary, but I haven't gone back and checked properly yet

def set_all_h_explicit(mol):

    for atom in mol.GetAtoms():

        atom.SetNumExplicitHs(atom.GetTotalNumHs())

        atom.SetNoImplicit(True)

# set_all_h_explicit

######

def apply(mol, balance_quat_surplus=False):

    mol = copy.deepcopy(mol)

    set_all_h_explicit(mol)

    pos  = [x[0] for x in mol.GetSubstructMatches(pos_pat)]
    quat = [x[0] for x in mol.GetSubstructMatches(quat_pat)]
    neg  = [x[0] for x in mol.GetSubstructMatches(neg_pat)]
    acid = [x[0] for x in mol.GetSubstructMatches(acid_pat)]

    logger.debug("{n_pos} positive/H, {n_quat} positive/quat and {n_neg} negative (of which {n_acid} are acid) charges identified".format(n_pos=len(pos), n_quat=len(quat), n_neg=len(neg), n_acid=len(acid)))

    h_added = 0

    # Negative charges...

    if quat:

        neg_surplus = len(neg) - len(quat)  # i.e. 'surplus' negative charges

        if neg_surplus > 0 and acid:

            logger.warn("zwitterion with more negative charges than quaternary positive centres detected")

            while neg_surplus > 0 and acid:

                atom = mol.GetAtomWithIdx(acid.pop(0))

                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                atom.SetFormalCharge(atom.GetFormalCharge() + 1)

                h_added += 1

                neg_surplus -= 1

        if balance_quat_surplus:

            quat_surplus = len(quat) - len(neg)

            acid_h = [x[0] for x in mol.GetSubstructMatches(acid_h_pat)]

            if quat_surplus > 0 and acid_h:

                logger.warn("Surplus of quat positive charges but with uncharged acids detected")

                while quat_surplus > 0 and acid_h:

                    atom = mol.GetAtomWithIdx(acid_h.pop(0))

                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() - 1)

                    h_added -= 1

                    quat_surplus -= 1

    else:

        for atom in [mol.GetAtomWithIdx(x) for x in neg]:

            while atom.GetFormalCharge() < 0:

                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)

                atom.SetFormalCharge(atom.GetFormalCharge() + 1)

                h_added += 1

    # Positive charges...

    for atom in [mol.GetAtomWithIdx(x) for x in pos]:

        while atom.GetFormalCharge() > 0 and atom.GetNumExplicitHs() > 0:

            atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)

            atom.SetFormalCharge(atom.GetFormalCharge() - 1)

            h_added -= 1

    # Done...

    try:

        sanity_check(mol)

    except StandardiseException as err:

        logger.debug("Molecule failed sanity check")

        raise

    logger.debug("Overall H balance: {sign}{n}; formal charge: {chg}".format(sign="+" if h_added > 0 else "", n=h_added, chg=formal_charge(mol)))

    return mol

# apply

####################################################################################################
# End
####################################################################################################
