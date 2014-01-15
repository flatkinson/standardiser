#! /usr/bin/env python

# Utility script to extract only those records where rules have been applied from SMILES output of standardiser.py

from __future__ import print_function

import sys

with open("standardised.smi", "r") as infile:
    with open("rules_applied.smi", "w") as outfile:

        for line in [x.strip() for x in infile]:

            if len(line.split()) == 4:

                print(line, file=outfile)
