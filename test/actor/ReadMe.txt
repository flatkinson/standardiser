# Structures from the EPA ACToR database
# http://actor.epa.gov/actor/faces/ACToRHome.jsp

PATH=../../bin:$PATH

./unique_smiles.py

wc -l actor.smi
  418667 actor.smi

PYTHONPATH=../.. standardiser.py -r actor.smi 1> standardiser.log 2>&1

rules_applied.py

# See notebook 'analysis.ipynb'
