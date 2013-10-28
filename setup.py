# Setup for the IPython Notebooks

import os, sys
import logging
from IPython.display import HTML
from jinja2 import Template
import pandas as pd
from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import AllChem, PandasTools, Draw
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.molSize = (450, 200)
