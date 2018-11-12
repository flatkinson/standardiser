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

import os
import sys
import argparse

import re
import csv
import json
from collections import Counter

from standardiser import standardise, SDF, make_logger
from standardiser.utils import errors

####################################################################################################

def main():

    ########################################################################
    # 
    # Program Parameters...
    # 

    script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]

    ######

    # Options, arguments and logging...

    argparser = argparse.ArgumentParser(description="Standardise compounds")

    argparser.add_argument("-V", "--verbose", action="store_true", help="enable verbose logger")
    argparser.add_argument("-r", "--output_rules_applied", action="store_true", help="enable output of rules applied")

    argparser.add_argument("-i", dest="infile", help="Input file (SDF or SMILES)")
    argparser.add_argument("-o", dest="outfile", help="Output file")

    config = argparser.parse_args()

    logger = make_logger.run(__name__)

    ######

    # Initialisation...

    rule_names = ["{:02d} {}".format(x['n'], x['name']) for x in standardise.rules.rule_set]

    counts = Counter({x: 0 for x in list(errors.keys()) + ['read', 'standardised']}) 

    input_type = os.path.splitext(config.infile)[1] # sdf or smi
    outfile_basename = os.path.splitext(config.infile)[0]
    outfile_ext      = os.path.splitext(config.infile)[1]

    ######

    logger.info("Input type = '{in_type}'".format(in_type=input_type))

    if input_type == ".sdf": # Read/write SDF...

        infile = SDF.readFile(open(config.infile))

        outfile = open(config.outfile, "w")
        errfile = open(outfile_basename + "_errors." + outfile_ext, "w")

        for original in infile:

            counts["read"] += 1

            logger.info(">>> Starting mol '{name}'...".format(name=original.name))

            ok = True

            try:

                if config.output_rules_applied:

                    rules_applied = []

                    parent = standardise.run(original.molblock, output_rules_applied=rules_applied)

                else:

                    parent = standardise.run(original.molblock)

            except standardise.StandardiseException as err:

                logger.warn(">>> {error} for '{name}'".format(error=errors[err.name], name=original.name))

                counts[err.name] += 1

                errfile.write("{mol}>  <n>\n{nread}\n\n>  <error>\n{error}\n\n$$$$\n".format(mol=original.molblock, nread=counts["read"], error=errors[err.name]))

                ok = False

            if ok:

                logger.info("Mol '{name}' OK".format(name=original.name))

                counts["standardised"] += 1

                parent = re.sub(r'^\w*\n', original.name + '\n', parent)

                if config.output_rules_applied:

                    rules_applied = ';'.join(rule_names[x-1] for x in rules_applied) if rules_applied else ''

                    outfile.write("{mol}>  <n>\n{nread}\n\n<rules_applied>\n{rules}\n\n$$$$\n".format(mol=parent, nread=counts["read"], rules=rules_applied))

                else:

                    outfile.write("{mol}>  <n>\n{nread}\n\n$$$$\n".format(mol=parent, nread=counts["read"]))

            if counts["read"] % 100 == 0: logger.info("...done: {read} read, {standardised} OK...".format(**counts))

    else: # Read/write (tab-seperated) SMILES + name...

        infile = csv.reader(open(config.infile), delimiter="\t")
        outfile = csv.writer(open(config.outfile, "w"), delimiter="\t")
        errfile_name = outfile_basename + "_errors." + outfile_ext
        errfile = csv.writer(open(errfile_name, "w"), delimiter="\t")

        for original in infile:

            counts["read"] += 1

            smiles, name = original

            logger.info(">>> Starting mol '{name}'...".format(name=name))

            ok = True

            try:

                if config.output_rules_applied:

                    rules_applied = []

                    parent = standardise.run(smiles, output_rules_applied=rules_applied)

                else:

                    parent = standardise.run(smiles)

            except standardise.StandardiseException as err:

                logger.warn(">>> {error} for mol '{name}'".format(error=errors[err.name], name=name))

                counts[err.name] += 1

                errfile.writerow(original + [err.name])

                ok = False

            if ok:

                logger.info("Mol '{name}' OK".format(name=name))

                counts["standardised"] += 1

                if config.output_rules_applied:

                    rules_applied = ';'.join(rule_names[x-1] for x in rules_applied) if rules_applied else ''

                    outfile.writerow([parent, name, smiles, rules_applied])

                else:

                    outfile.writerow([parent, name])

            if counts["read"] % 100 == 0: logger.info("...done: {read} read, {standardised} OK...".format(**counts))

    logger.info("Finished: {read} read, {standardised} OK in total.".format(**counts))

    logger.info("Counts: " + json.dumps(counts, indent=4))

# main

####################################################################################################
    
if __name__ == '__main__':
    main()    

####################################################################################################
# End
####################################################################################################
