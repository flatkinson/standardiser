#! /usr/bin/env python

import csv
import random

sample_size = 10000

with open("chembl_parents.smi", "r") as csv_file:

    reader = csv.reader(csv_file, delimiter="\t")

    smiles = list(reader)

subset = [smiles[x] for x in random.sample(xrange(len(smiles)), sample_size)]

with open("subset.smi", "w") as csv_file:

    writer = csv.writer(csv_file, delimiter="\t")

    for record in subset:

        writer.writerow(record)
