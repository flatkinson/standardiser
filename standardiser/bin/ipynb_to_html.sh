#! /bin/bash

mkdir html

\ls -1 *.ipynb | while read ipynb; do stem=`basename ${ipynb} \.ipynb`; jupyter nbconvert --to=html --output html/${stem}.html ${ipynb}; done

cd html/
perl -i -p -e 's#(<a href=".+?)\.ipynb#$1.html#g' *.html
cp ../standardiser.pdf .
cd ..
