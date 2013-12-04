# Check bahaviour of standardiser by running on ChEMBL parent structures

./chembl_parents.py

./subset.py

rm -f not_built.smi multi_component.smi no_non_salt.smi standardised.smi

../standardiser2.py subset.smi 1> standardiser.log 2>&1

# $ wc -l not_built.smi multi_component.smi standardised.smi no_non_salt.smi
#      0 not_built.smi
#      0 multi_component.smi
#   9994 standardised.smi
#      6 no_non_salt.smi
#  10000 total

# Start ipython notebook server and look at 'chembl_parents.ipynb'
