#! /usr/bin/env python -u

# subset a fraction of the lines in a file

import sys, os
import random

chunk_size = 100

input_filename = sys.argv[1]
sample_fraction = float(sys.argv[2])

sample_size = int(chunk_size*sample_fraction)

stem, ext = os.path.splitext(os.path.basename(input_filename))

with open(input_filename, "r") as input_file:

    with open("{}_{:03d}{}".format(stem, int(round(sample_fraction*100)), ext), "w") as output_file:

        records = []

        for n, record in enumerate(input_file, 1):

            records.append(record)

            if n % chunk_size == 0:

                output_file.writelines(records[x] for x in random.sample(list(range(chunk_size)), sample_size))

                records = []

        output_file.writelines(records[x] for x in random.sample(list(range(len(records))), int(round(len(records)*sample_fraction))))
