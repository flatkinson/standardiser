import sys, os, logging

### logging.basicConfig(level=logging.INFO, format="[%(asctime)s %(levelname)-8s %(name)s] %(message)s", datefmt="%Y/%b/%d %H:%M:%S")
logging.basicConfig(level=logging.INFO, format="[%(asctime)s %(levelname)-8s] %(message)s", datefmt="%Y/%b/%d %H:%M:%S")

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Geometry import rdGeometry
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.molSize = (450, 200)

sys.path.append("/Users/francis/projects/standardiser")
