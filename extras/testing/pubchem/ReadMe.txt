# Structures from PubChem
# ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/

wc -l pubchem.smi
 47332186 pubchem.smi

subset2.py pubchem.smi 0.02

wc -l pubchem_002.smi
  946643 pubchem_002.smi

PYTHONPATH=../.. ../../bin/standardiser.py -r pubchem_002.smi 1> standardiser.log 2>&1

rules_applied.py

# See notebook 'analysis.ipynb' 
