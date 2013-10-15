#! /bin/bash

rm -f *.slides.html standardiser.html

for nb in *.ipynb
do
	ipython nbconvert --to slides ${nb}
done

(
	echo "<h3>Demonstration of standardiser modules<h3>"

	for page in *.slides.html
	do
		name=`basename ${page} '.slides.html'`

		echo "<li><a target="_blank" href="http://127.0.0.1:8000/${page}">${name}</a></li>"
	done

	echo "</li>"

) > standardiser.html

python -m SimpleHTTPServer 1> SimpleHTTPServer.log 2>&1 &

# open -a /Applications/Safari.app standardiser.html
open -a /Applications/Google\ Chrome.app standardiser.html
