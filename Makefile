build: doc
	R -e "devtools::build()"

doc:
	R -e "devtools::document()"
	R -e "devtools::build_vignettes()"

check: 
	R -e "devtools::check()"

install: 
	R CMD INSTALL .

first.install:
	mv DESCRIPTION DESCRIPTION.save
	sh dev/install_ubuntu.sh
	Rscript dev/create_description.R
	R -e "devtools::install_local('.',dependencies=TRUE, force=TRUE)"

VERSION = "v3.5.0-beta"

tags:
	git tag -f -a $(VERSION) -m 'version $(VERSION)'
	git push && git push --tags

clean:
	rm `find . -name '.#*' -o -name '#*' -o -name '*~'`

# filter by file name (eg. ensembl) will run tests inside file test-pgx-ensembl.R
filter=
test:
	R -e "devtools::test(filter='$(filter)')"

FORCE: ;

docker.os: FORCE
	docker build --progress plain \
	  -f dev/Dockerfile.os -t playbase-os . \
	2>&1 | tee docker-os.log

docker: FORCE
	docker build --progress plain \
	  -f dev/Dockerfile -t bigomics/playbase . \
	2>&1 | tee docker.log
