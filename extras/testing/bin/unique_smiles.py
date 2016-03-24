#! /usr/bin/env python

import sys, os
import csv, collections

unique = collections.defaultdict(list)

# Read SMILES + name...

with open("actor.csv") as infile:

    for record in csv.DictReader(infile, delimiter=","):

        unique[record["smiles"]].append(record["cas_number"])

# Output unique SMILES + list of names...

with open("actor.smi", "w") as outfile:

    writer = csv.writer(outfile, "w", delimiter="\t")

    for smiles, names in unique.items():

        writer.writerow([smiles, ":".join(names)])
