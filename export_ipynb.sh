#! /bin/bash

file=index.html

echo "<ul>" > ${file}

for ipynb in *.ipynb
do

	name=`basename ${ipynb} '.ipynb'`

	ipython nbconvert ${ipynb}

	echo "<li><a target="_blank" href="${name}.html">${name}</a></li>"

done >> ${file}

echo "</ul>" >> ${file}

# open -a /Applications/Google\ Chrome.app ${file}
