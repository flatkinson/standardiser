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
This module provides a lightweight utility to read molecule blocks from an SD file. This is so the block can be
written to a file for later inspection if RDKit cannot build a molecule from it. The readers in RDKit return
nothing if they cannot build the molecule, so the molblock cannot easily be saved.
"""

####################################################################################################

import re

####################################################################################################

def readFile(filename):

    lines = []

    for line in open(filename) if type(filename) == str else filename:

        lines.append(line)

        if line == "$$$$\n":

            molblock = "".join(lines)

            lines = []

            yield Molfile(molblock)

class Molfile(dict):

    n_mols = 0

    def __init__(self, molblock):

        self.__class__.n_mols += 1

        self.original = molblock

        self.molblock, data = re.search("\A(.*\nM\s+END\s*\n)(.*)\$\$\$\$", molblock, re.DOTALL).groups()

        data = dict(re.findall("^>\s+<(.*?)>(?:\s*\(\d+\))?\s*\n(.*?)\s*\n\n", data, re.DOTALL | re.MULTILINE))

        self.name = re.match('^(\w+)?\n', molblock).group(1) or data.get("Name") or data.get("name") or data.get("molregno") or "mol_{n:04d}".format(n=self.__class__.n_mols)

        self.__dict__.update(data)

    def __getitem__(self, key):

        return self.__dict__[key]

    def write(self, name="junk"):

        open(name + ".mol", "w").write(self.original)

####################################################################################################
# End
####################################################################################################
