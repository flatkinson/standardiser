# Copy HTML to wwwdev...

for dir in docs issues test; do ( cd $dir; ../export_ipynb.sh; ) done

ssh ebi-001 "cd /homes/francis/wwwdev && mkdir standardiser && cd standardiser && mkdir docs test"

rsync -avz index.html standardiser.pptx ebi-001:/homes/francis/wwwdev/standardiser

rsync -avz docs/*.html ebi-001:/homes/francis/wwwdev/standardiser/docs
rsync -avz test/*.html ebi-001:/homes/francis/wwwdev/standardiser/test

# Make tarball...

rm -f docs/*.html test/*.html

cd ..

tar -czvf standardiser.tgz \
	--exclude='archive.sh' \
	--exclude='problems' \
	--exclude='OLD' \
	--exclude='*.pyc' \
	--exclude='*.log' \
	--exclude='.git*' \
	--exclude='.ipynb_checkpoints' \
	--exclude='*.swp' \
	standardiser/
