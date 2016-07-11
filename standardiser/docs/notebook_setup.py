# Common setup for the IPython Notebooks

from __future__  import print_function, division, absolute_import
import six

import warnings

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
