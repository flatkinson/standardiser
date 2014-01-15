Investigation of some problem structures from Ulf Norinder

Ten structures contain 'radical' groups (i.e. the molblock includes a 'RAD' tag). RDKit will not build these when reading SDF format (although it accespts them as SMILES).
These structures are in not_built.sdf, and are depicted in not_built_from_sdf.png (via MarvinView).

Of these, nine are simply metal counterions that appear to have been 'neutralised' in a way that RDKit does not understand.

One, 377718, appears to contain a guanidiumium moiety that may have been improperly neutralised.

For the investigation in the IPython Notebook 'problems.ipynb', the radical flags were removed from problems.sdf (the original SDF file is 'problems.sdf.01').
