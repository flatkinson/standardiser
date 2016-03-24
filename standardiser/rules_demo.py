import os
import re
from jinja2 import Template
from tempfile import TemporaryFile
from base64 import b64encode

from rdkit import Chem
from rdkit.Chem import Draw

from . import rules

####################################################################################################

# Dir for Jinja2 templates...

templates = os.path.join(os.path.dirname(__file__), 'templates')

####################################################################################################

def rules_table():

    def f(x):

        x = x.copy()
            
        x['tag'] = re.sub("[ ,]", "_", re.sub("[()>.]", "", x["name"]))
            
        return x
                    
    rule_set = [f(x) for x in rules.rule_set]

    template = Template(open(os.path.join(templates, 'rules_table.html')).read())

    return template.render(rule_set=rule_set)

######

def b64_img(mol, match):
    
    with TemporaryFile() as fh:

        Draw.MolToImage(mol, highlightAtoms=match).save(fh, format='png')

        fh.seek(0)

        b64_img = b64encode(fh.read()).decode('utf-8')
        
    return b64_img

######

def show_change(smiles):

    mol = Chem.MolFromSmiles(smiles)

    old_mol, old_match, new_mol, new_match = rules.demo(mol)

    template = Template(open(os.path.join(templates,  'show_rule.html')).read())

    if not new_mol: return template.render(old_img=b64_img(mol, None), new_img=None)

    return template.render(old_img=b64_img(old_mol, old_match), new_img=b64_img(new_mol, new_match))

####################################################################################################
# End
####################################################################################################
