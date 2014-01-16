# Put setup for the IPython Notebook here to reduce clutter at top of page

from __future__  import print_function, division

import os, sys, logging
import re
import random

from jinja2 import Template

import pandas as pd
import numpy as np

import psycopg2, psycopg2.extras

from IPython.display import HTML

from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import AllChem, PandasTools, Draw
from rdkit.Chem.Draw import IPythonConsole

PandasTools.RenderImagesInAllDataFrames()
IPythonConsole.molSize = (450, 200)

# At top of notebook...
# %run setup.py
# 
# reload(logging) # Seems to be required to get logging from modules working
# 
# logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s %(levelname)-8s] %(message)s", datefmt="%Y/%b/%d %H:%M:%S")
