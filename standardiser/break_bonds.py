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
Break bonds to Group I and II metals
"""

####################################################################################################

from . import make_logger
logger = make_logger.run(__name__)

from collections import defaultdict

from rdkit import Chem

from .utils import StandardiseException, sanity_check

####################################################################################################

# Module configuration...

bonds_to_break = [("[Li,Na,K,Mg,Ca]-[#7,#8,#16]", 1), ("[Mg,Ca]=N", 2)]  # SMARTS pattern, charge increment

####################################################################################################

# Module initialization...

bonds_to_break = [(Chem.MolFromSmarts(smarts), charge_incr) for smarts, charge_incr in bonds_to_break]

####################################################################################################

def run(mol):

    charge_added = defaultdict(int)

    n_broken = 0

    for pattern, charge_incr in bonds_to_break:

        idx_pairs = mol.GetSubstructMatches(pattern)

        ed_mol = Chem.EditableMol(mol)

        for i, j in idx_pairs:

            ed_mol.RemoveBond(i, j)

            n_broken += 1

        mol = ed_mol.GetMol()

        for i, j in idx_pairs:

            atom = mol.GetAtomWithIdx(i)
            atom.SetFormalCharge(atom.GetFormalCharge() + charge_incr)

            charge_added[i] += charge_incr

            atom = mol.GetAtomWithIdx(j)
            atom.SetFormalCharge(atom.GetFormalCharge() - charge_incr)

            charge_added[j] -= charge_incr

    try:

        sanity_check(mol)

    except StandardiseException as err:

        logger.debug("Molecule failed sanity check")

        raise

    logger.debug("Broke {n} bonds to Group I and II metals".format(n=n_broken))

    return mol

# run

####################################################################################################
# End
####################################################################################################
