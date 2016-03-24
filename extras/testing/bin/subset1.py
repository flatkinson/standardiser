#! /usr/bin/env python

# subset a fixed number of lines in a file

import sys, os
import random

input_filename = sys.argv[1]
sample_size = int(sys.argv[2])

stem, ext = os.path.splitext(os.path.basename(input_filename))

with open(input_filename, "r") as input_file:

    records = input_file.readlines()

subset = [records[x] for x in random.sample(list(range(len(records))), sample_size)]

with open("{}_{}{}".format(stem, sample_size, ext), "w") as output_file:

    output_file.writelines(subset)
