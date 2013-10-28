# Check bahaviour of standardiser by running on ChEMBL parent structures

./chembl_parents.py

# $ wc -l chembl_parents.smi
#  10000 chembl_parents.smi

rm -f not_built.smi multi_component.smi no_non_salt.smi standardised.smi

../standardiser.py chembl_parents.smi

# $ wc -l not_built.smi multi_component.smi standardised.smi no_non_salt.smi
#      0 not_built.smi
#      0 multi_component.smi
#   9994 standardised.smi
#      6 no_non_salt.smi
#  10000 total

# ipython notebook
