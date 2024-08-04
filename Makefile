build: doc
	R -e "devtools::build()"

doc:
	R -e "devtools::document()"
	R -e "devtools::build_vignettes()"

check: 
	R -e "devtools::check()"

update: 
	R -e "source('dev/rspm.R');BiocManager::install(ask=FALSE)"

install: 
	R CMD INSTALL .

installx: 
	Rscript dev/install_playbase.R

full.install:
	sudo sh dev/install_ubuntu.sh
	Rscript dev/create_description.R
	Rscript dev/install_playground.R

VERSION="v3.5.0-beta"

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

docker.os: 
	docker build --no-cache \
	  -f dev/Dockerfile.os -t playbase-os . \
	  2>&1 | tee docker-os.log

docker: 
	docker build --no-cache \
	  -f dev/Dockerfile -t bigomics/playbase . \
	  2>&1 | tee docker.log
dockerx:
	cat dev/Dockerfile.os dev/Dockerfile > dev/Dockerfilex
	docker build --no-cache \
	  -f dev/Dockerfilex -t bigomics/playbase:x . \
	  2>&1 | tee dockerx.log
