import os, sys
import logging

import pandas as pd

from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import AllChem, PandasTools, Draw
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.molSize = (450, 200)
from IPython.display import HTML

sys.path.append("/Users/francis/projects/standardiser")
