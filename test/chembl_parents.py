#! /usr/bin/env python

# Get parent structures from CHEMBL as RDKit canonical SMILES

from __future__ import print_function, division

import cx_Oracle
from rdkit import Chem

maxrows = 10000

conn_str = 'chembl_17_app/chembl_17_app@chempro'

sql = """
select
      b.canonical_smiles
    , c.chembl_id
from
      (select distinct parent_molregno as molregno from molecule_hierarchy) a
    , compound_structures b
    , chembl_id_lookup c
    , compound_properties d
where
    a.molregno = b.molregno
and a.molregno = c.entity_id and c.entity_type = 'COMPOUND'
and a.molregno = d.molregno
and b.canonical_smiles is not null
and d.mw_freebase >= 100 and d.mw_freebase <= 400
and rownum <= :maxrows
"""

conn = cx_Oracle.connect(conn_str)

cursor = conn.cursor()

cursor.execute(sql, maxrows=maxrows)

output = open("chembl_parents.smi", "w")
errors = open("chembl_parents.err", "w")

for row in cursor:

    smiles, name = row

    mol = Chem.MolFromSmiles(smiles)

    if not mol:

        print("\t".join([smiles, name]), file=errors)

    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

    print("\t".join([smiles, name]), file=output)

output.close()
