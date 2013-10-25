#! /usr/bin/env python

from __future__ import print_function
import os, sys, argparse, logging

import copy

#@ from collections import Counter

from rdkit import Chem

from standardise import standardise

########################################################################
# 
# Program Parameters...
# 

script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]

######

# Options and arguments...

argparser = argparse.ArgumentParser(description="Standardize compounds")

argparser.add_argument("-V", "--verbose", action="store_true", help="enable verbose logging")

argparser.add_argument("infile", help="Input file (SDF or SMILES)")

config = argparser.parse_args()

######

# Initialisation...

logging.basicConfig(level=logging.DEBUG if config.verbose else logging.INFO, format="[%(asctime)s %(levelname)-8s] %(message)s", datefmt="%Y/%b/%d %H:%M:%S")

########################################################################

# Initialise...

error_names = standardise.errors.keys()

#@ counts = Counter({name: 0 for name in error_names + ["read", "standardised"]})
counts = dict((name, 0) for name in error_names + ["read", "standardised"]) # python2.6-compatible

# Input type...

input_type = os.path.splitext(config.infile)[1] # sdf or smi

logging.info("Input type = '{type}'".format(type=input_type))

if input_type == ".sdf":

    # Read/write sdf...

#@    outfile = {name: open(name + ".sdf", "w") for name in error_names + ["standardised"]}
    outfile = dict((name, open(name + ".sdf", "w")) for name in error_names + ["standardised"]) # python2.6-compatible

    for sdf in standardise.SDF.readFile(open(config.infile)):

        counts["read"] += 1

        logging.info(">>> Starting mol '{name}'...".format(name=sdf.name))

        try:

            parent = standardise.apply(sdf.molblock)

        except standardise.StandardiseException as e:

            logging.warn(">>> {err} for '{name}'".format(err=standardise.errors[e.name], name=sdf.name))

            counts[e.name] += 1

            outfile[e.name].write(sdf.original)

            continue

        logging.info("Mol '{name}' OK".format(name=sdf.name))

        counts["standardised"] += 1

        outfile["standardised"].write("{molblock}>  <n>\n{n_read}\n\n$$$$\n".format(molblock=parent, n_read=counts["read"]))

else: # input_type == '.smi'

    # Read/write smi...

    #@ outfile = {name: open(name + ".smi", "w") for name in error_names + ["standardised"]}
    outfile = dict((name, open(name + ".smi", "w")) for name in error_names + ["standardised"])  # python2.6-compatible

    for original in open(config.infile):

        smiles, name = original.split()

        counts["read"] += 1

        logging.info(">>> Starting mol '{name}'...".format(name=name))

        try:

            parent = standardise.apply(smiles)

        except standardise.StandardiseException as e:

            logging.warn(">>> {err} for '{name}'".format(err=standardise.errors[e.name], name=name))

            counts[e.name] += 1

            outfile[e.name].write(original)

            continue

        logging.info("Mol '{name}' OK".format(name=name))

        counts["standardised"] += 1

        outfile["standardised"].write("{smiles}\t{name}\n".format(smiles=parent, name=name))

logging.info("Finished: {read} read, {not_built} not built, {no_non_salt} no non-salt, {multi_component} multi-component, {standardised} standardised".format(**counts))

########################################################################
# End
########################################################################
