rm -f docs/*.html issues/*.html

cd ..

tar -czvf standardiser.tgz \
	--exclude 'problems' \
	--exclude='OLD' \
	--exclude '*.pyc' \
	--exclude='*.log' \
	--exclude '.git*' \
	--exclude='.ipynb_checkpoints' \
	--exclude='*.swp' \
	standardiser/
