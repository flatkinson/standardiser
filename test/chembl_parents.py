#! /usr/bin/env python

# Get parent structures from CHEMBL

from __future__ import print_function, division

import psycopg2
import csv

db = {
        "host":   "10.7.248.39"
      , "port":   5432
      , "user":   "chembl"
      , "dbname": "chembl_17"
}

sql = """
select
      b.canonical_smiles as smiles
    , c.chembl_id
from
      (select distinct parent_molregno as molregno from molecule_hierarchy) a
    , compound_structures b
    , chembl_id_lookup c
    , compound_properties d
    , mols_rdkit e
where
    a.molregno = b.molregno
and a.molregno = c.entity_id and c.entity_type = 'COMPOUND'
and a.molregno = d.molregno
and a.molregno = e.molregno
and b.canonical_smiles is not null
and e.m is not null
and d.heavy_atoms >= 10 and d.heavy_atoms <= 40
and d.mw_freebase <= 500
"""

with psycopg2.connect(**db) as conn:

    cursor = conn.cursor()

    cursor.execute(sql)

    with open("chembl_parents2.smi", "w") as tsv_file:

        writer = csv.writer(tsv_file, delimiter="\t")

        for record in cursor:

            writer.writerow(record)
