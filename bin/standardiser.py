#! /usr/bin/env python -u

# Simple driver program for standardise module
# 
# TODO: get rid of duplicated code

########################################################################

import os, sys, argparse, logging

import csv, json

from rdkit import Chem

from standardise import standardise, SDF
from standardise.utils import errors

########################################################################
# 
# Program Parameters...
# 

script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]

######

# Options and arguments...

argparser = argparse.ArgumentParser(description="Standardise compounds")

argparser.add_argument("-V", "--verbose", action="store_true", help="enable verbose logging")
argparser.add_argument("-r", "--output_rules_applied", action="store_true", help="enable output of rules applied")

argparser.add_argument("infile", help="Input file (SDF or SMILES)")

config = argparser.parse_args()

######

# Initialisation...

logging.basicConfig(level=logging.DEBUG if config.verbose else logging.INFO, format="[%(asctime)s %(levelname)-8s] %(message)s", datefmt="%Y/%b/%d %H:%M:%S")

########################################################################

# Initialise...

error_names = errors.keys()

counts = dict((name, 0) for name in error_names + ["read", "standardised"]) # python2.6-compatible

# Input type...

input_type = os.path.splitext(config.infile)[1] # sdf or smi

logging.info("Input type = '{}'".format(input_type))

if input_type == ".sdf": # Read/write SDF...

    infile = SDF.readFile(open(config.infile))

    outfile = open("standardised.sdf", "w")
    errfile = open("errors.sdf", "w")

    for original in infile:

        counts["read"] += 1

        logging.info(">>> Starting mol '{}'...".format(original.name))

        ok = True

        try:

            if config.output_rules_applied:

                rules_applied = []

                parent = standardise.apply(original.molblock, output_rules_applied=rules_applied)

            else:

                parent = standardise.apply(original.molblock)

        except standardise.StandardiseException as err:

            logging.warn(">>> {} for '{}'".format(errors[err.name], original.name))

            counts[err.name] += 1

            errfile.write("{}>  <n>\n{}\n\n<error>\n{}\n\n$$$$\n".format(original.molblock, counts["read"], errors[err.name]))

            ok = False

        if ok:

            logging.info("Mol '{}' OK".format(original.name))

            counts["standardised"] += 1

            if config.output_rules_applied:

                outfile.write("{}>  <n>\n{}\n\n<rules_applied>\n{}\n\n$$$$\n".format(parent, counts["read"], ','.join(str(x) for x in rules_applied)))

            else:

                outfile.write("{}>  <n>\n{}\n\n$$$$\n".format(parent, counts["read"]))

        if counts["read"] % 100 == 0: logging.info("...done: {read} read, {standardised} OK...".format(**counts))

else: # Read/write (tab-seperated) SMILES + name...

    infile = csv.reader(open(config.infile), delimiter="\t")

    outfile = csv.writer(open("standardised.smi", "w"), delimiter="\t")
    errfile = csv.writer(open("error.smi", "w"), delimiter="\t")

    for original in infile:

        counts["read"] += 1

        smiles, name = original

        logging.info(">>> Starting mol '{}'...".format(name))

        ok = True

        try:

            if config.output_rules_applied:

                rules_applied = []

                parent = standardise.apply(smiles, output_rules_applied=rules_applied)

            else:

                parent = standardise.apply(smiles)

        except standardise.StandardiseException as err:

            logging.warn(">>> {} for mol '{}'".format(errors[err.name], name))

            counts[err.name] += 1

            errfile.writerow(original + [err.name])

            ok = False

        if ok:

            logging.info("Mol '{}' OK".format(name))

            counts["standardised"] += 1

            if config.output_rules_applied:

                outfile.writerow([parent, name, smiles, ','.join(str(x) for x in rules_applied)])

            else:

                outfile.writerow([parent, name])

        if counts["read"] % 100 == 0: logging.info("...done: {read} read, {standardised} OK...".format(**counts))

logging.info("Finished: {read} read, {standardised} OK in total.".format(**counts))

logging.info("Counts: " + json.dumps(counts, indent=4))

########################################################################
# End
########################################################################
