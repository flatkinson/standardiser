# Common setup for the IPython Notebooks

from __future__  import print_function, division, absolute_import
import six

import warnings
import logging

# NB These steps together seem to be required to get logging from modules working, but I don't know why yet...
if six.PY3: from importlib import reload
reload(logging) 
logging.basicConfig(level=logging.INFO, format="[%(asctime)s %(levelname)-8s] %(message)s", datefmt="%Y/%b/%d %H:%M:%S")

import sys
import re
import random
 
from ipywidgets import HTML

import pandas as pd
 
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools, AllChem
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.molSize = (450, 200)
PandasTools.RenderImagesInAllDataFrames()
