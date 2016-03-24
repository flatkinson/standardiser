#! /usr/bin/env python
####################################################################################################
# 
# Copyright [2014] EMBL - European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in
# compliance with the License.  You may obtain a copy of
# the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software distributed under the
# License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
# either express or implied. See the License for the specific language governing permissions
# and limitations under the License.
# 
####################################################################################################

# Simple driver program for standardise module
# 
# TODO: get rid of duplicated code

####################################################################################################

from __future__ import print_function, division, absolute_import
import six

import os
import sys
import argparse
import logging

import csv
import json
from collections import Counter

from rdkit import Chem

from standardiser import standardise, SDF
from standardiser.utils import errors

####################################################################################################

def main():

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

    counts = Counter(list(errors.keys()) + ['read', 'standardised']) 

    # Input type...

    input_type = os.path.splitext(config.infile)[1] # sdf or smi

    logging.info("Input type = '{in_type}'".format(in_type=input_type))

    if input_type == ".sdf": # Read/write SDF...

        infile = SDF.readFile(open(config.infile))

        outfile = open("standardised.sdf", "w")
        errfile = open("errors.sdf", "w")

        for original in infile:

            counts["read"] += 1

            logging.info(">>> Starting mol '{name}'...".format(name=original.name))

            ok = True

            try:

                if config.output_rules_applied:

                    rules_applied = []

                    parent = standardise.apply(original.molblock, output_rules_applied=rules_applied)

                else:

                    parent = standardise.apply(original.molblock)

            except standardise.StandardiseException as err:

                logging.warn(">>> {error} for '{name}'".format(error=errors[err.name], name=original.name))

                counts[err.name] += 1

                errfile.write("{mol}>  <n>\n{nread}\n\n<error>\n{error}\n\n$$$$\n".format(mol=original.molblock, nread=counts["read"], error=errors[err.name]))

                ok = False

            if ok:

                logging.info("Mol '{name}' OK".format(name=original.name))

                counts["standardised"] += 1

                if config.output_rules_applied:

                    outfile.write("{mol}>  <n>\n{nread}\n\n<rules_applied>\n{rules}\n\n$$$$\n".format(mol=parent, nread=counts["read"], rules=','.join(str(x) for x in rules_applied)))

                else:

                    outfile.write("{mol}>  <n>\n{nread}\n\n$$$$\n".format(mol=parent, nread=counts["read"]))

            if counts["read"] % 100 == 0: logging.info("...done: {read} read, {standardised} OK...".format(**counts))

    else: # Read/write (tab-seperated) SMILES + name...

        infile = csv.reader(open(config.infile), delimiter="\t")

        outfile = csv.writer(open("standardised.smi", "w"), delimiter="\t")
        errfile = csv.writer(open("errors.smi", "w"), delimiter="\t")

        for original in infile:

            counts["read"] += 1

            smiles, name = original

            logging.info(">>> Starting mol '{name}'...".format(name=name))

            ok = True

            try:

                if config.output_rules_applied:

                    rules_applied = []

                    parent = standardise.apply(smiles, output_rules_applied=rules_applied)

                else:

                    parent = standardise.apply(smiles)

            except standardise.StandardiseException as err:

                logging.warn(">>> {error} for mol '{name}'".format(error=errors[err.name], name=name))

                counts[err.name] += 1

                errfile.writerow(original + [err.name])

                ok = False

            if ok:

                logging.info("Mol '{name}' OK".format(name=name))

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
    
if __name__ == '__main__':
    main()    

####################################################################################################
# End
####################################################################################################
