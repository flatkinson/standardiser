#! /bin/bash

rm -f *.html *.pdf

ipython nbconvert ../*.ipynb

perl -i -p -e 's/\.ipynb/.html/g' *.html

perl -i -p -e 's#(?<=src=")files/##g' *.html

ln -s ../standardiser.pdf .

ln -s 00_Introduction.html index.html
