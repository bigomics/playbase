build: doc
	R -e "devtools::build()"

doc:
	R -e "devtools::document()"
	R -e "devtools::build_vignettes()"

check: 
	R -e "devtools::check()"

install: 
	R CMD INSTALL .

install.dep:
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

docker: FORCE
	docker build -f dev/Dockerfile -t playbase .

docker.nc: FORCE
	docker build --no-cache -f dev/Dockerfile -t playbase .

FORCE: ;
