import urllib
import StringIO
from IPython.display import HTML
from jinja2 import Template

from rdkit import Chem
from rdkit.Chem import Draw

from standardise import rules

def img(mol, match):

    temp = StringIO.StringIO()

    Draw.MolToImage(mol, highlightAtoms=match).save(temp, format="PNG")

    img = temp.getvalue().encode('base64')

    temp.close()

    return img

# img_tag

def demo(smiles):

    mol = Chem.MolFromSmiles(smiles)

    old_mol, old_match, new_mol, new_match = rules.demo(mol)

    return HTML(Template(open("templates/old_vs_new.html").read()).render(old_img=img(old_mol, old_match), new_img=img(new_mol, new_match)))

# demo

def s2m(smiles):

    return Chem.MolFromSmiles(smiles)

# s2m
