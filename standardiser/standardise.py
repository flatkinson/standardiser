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
Apply standardisation procedure
"""

####################################################################################################

from . import make_logger
logger = make_logger.run(__name__)

from rdkit import Chem

from . import break_bonds, neutralise, rules, unsalt

from .utils import StandardiseException, sanity_check, timeout

####################################################################################################
#
# Module configuration...
#

####################################################################################################
#
# Module initialization...
#

####################################################################################################

def verbose(verbose=True):

    """
    Turn full debugging output on or off.
    """

    level = 10 if verbose else 20

    logger.setLevel(level)

    for module in break_bonds, unsalt, neutralise, rules: module.logger.setLevel(level)

####################################################################################################

# Unix signals such as SIGALRM are not unavailable on Windows, so the timeout facility cannot be used.
# Use of signals such as SIGALRM is also impossible whe running under mod_wsgi...
# 	https://github.com/GrahamDumpleton/mod_wsgi-docs/blob/master/configuration-directives/WSGIRestrictSignal.rst
# 	https://github.com/GrahamDumpleton/mod_wsgi-docs/blob/master/developer-guides/tips-and-tricks.rst
# Thus, the use of the global use of the timeout wrapper is disabled here and enabled selectively below.

### @timeout()
def run(input_mol, output_rules_applied=None, verbose=False): 

    # Get input molecule...

    if type(input_mol) == Chem.rdchem.Mol:

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

    try:

        sanity_check(mol)

    except StandardiseException as err:

        logger.debug("Molecule failed sanity check")

        raise

    ######

    # Get disconnected fragments...

    non_salt_frags = []

    mol = break_bonds.run(mol)

    for n, frag in enumerate(Chem.GetMolFrags(mol, asMols=True), 1):

        logger.debug("Starting fragment {n} '{smi}'...".format(n=n, smi=Chem.MolToSmiles(frag)))

        logger.debug("1) Check for non-organic elements...")

        if unsalt.is_nonorganic(frag): continue

        logger.debug("2) Attempting to neutralise (first pass)...")

        frag = neutralise.run(frag)

        logger.debug("3) Applying rules...")

        frag = rules.run(frag, output_rules_applied=output_rules_applied, verbose=verbose)

        logger.debug("4) Attempting to neutralise (second pass)...")

        frag = neutralise.run(frag)

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

# run

######

# Check for availability of timeout before enabling (see above for details)...

import platform

if platform.system() == 'Windows':

    logger.warning("Running under Windows: must disable use of timeout")

else:

    try:

        from mod_wsgi import version

        logger.warning("Running under mod_wsgi: must disable use of timeout")

    except:

        run = timeout()(run)

####################################################################################################
# End
####################################################################################################
