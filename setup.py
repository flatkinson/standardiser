import os, sys
import logging

from IPython.display import HTML

import pandas as pd

from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import AllChem, PandasTools, Draw
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.molSize = (450, 200)
