#! /bin/bash

file=issues.html

echo "<h2>Discussion of some issues around <code>standardiser</code> modules<h2>" > ${file}

echo "<ul>" >> ${file}

for ipynb in *.ipynb
do

	name=`basename ${ipynb} '.ipynb'`

	ipython nbconvert ${ipynb}

	echo "<li><a target="_blank" href="${name}.html">${name}</a></li>"

done >> ${file}

echo "</ul>" >> ${file}

open -a /Applications/Google\ Chrome.app ${file}
