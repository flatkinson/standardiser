# Setup for the IPython Notebooks

import os, sys
import logging
from IPython.display import HTML
from jinja2 import Template
import pandas
import pandas as pd
from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import AllChem, PandasTools, Draw
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.molSize = (450, 200)

# At top of notebook...
# %run setup.py
# 
# reload(logging) # Seems to be required to get logging from modules working
# 
# logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s %(levelname)-8s] %(message)s", datefmt="%Y/%b/%d %H:%M:%S")
