"""SD file  utilities"""

import re
from rdkit import Chem

def readFile(filename):

    lines = []

    for line in open(filename) if type(filename) == str else filename:

        lines.append(line)

        if line == "$$$$\n":

            molfile = "".join(lines)

            lines = []

            yield Molfile(molfile)

class Molfile(dict):

    n_mols = 0

    def __init__(self, molfile):

        self.__class__.n_mols += 1

        self.original = molfile

        self.molblock, data = re.search("\A(.*\nM\s+END\s*\n)(.*)\$\$\$\$", molfile, re.DOTALL).groups()

        data = dict(re.findall("^>\s+<(.*?)>(?:\s*\(\d+\))?\s*\n(.*?)\s*\n\n", data, re.DOTALL | re.MULTILINE))

        self.__dict__.update(data)

        self.name = data.get("Name") or data.get("name") or data.get("molregno") or "mol_{:04d}".format(self.__class__.n_mols)

    def __getitem__(self, key):

        return self.__dict__[key]

    def write(self, name="junk"):

        open(name + ".mol", "w").write(self.original)

