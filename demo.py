import logging

from rdkit import Chem

from standardise import standardise

for input in [ "Oc1ncccc1.O.O", "CC(=O)O[Na].O.Cl" ]:

	print("Input: {}".format(input))

	mol = Chem.MolFromSmiles(input)

	molblock = Chem.MolToMolBlock(mol)

	try:

		parent = standardise.apply(molblock)

	except standardise.StandardiseException as e:

		print("Standardisation failed: {}".format(standardise.errors[e.name]))

		continue

	output = Chem.MolToSmiles(Chem.MolFromMolBlock(parent))

	print("Output: {}\n".format(output))
