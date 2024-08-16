BRANCH:=`git rev-parse --abbrev-ref HEAD`  ## get active GIT branch
BRANCH:=$(strip $(BRANCH))
TAG=$(BRANCH)
VERSION="v3.5.0-beta"
ARG=

build: doc
	R -e "devtools::build()"

doc:
	R -e "devtools::document()"
	R -e "devtools::build_vignettes()"

check: 
	R -e "devtools::check()"

depend: 
	Rscript dev/install_dependencies.R 

install: depend
	R CMD INSTALL .

install.rcmd: 
	Rscript dev/install_playbase.R 'rcmd'

install.local: 
	Rscript dev/install_playbase.R 'local'

install.full:
	sudo sh dev/install_ubuntu.sh
	Rscript dev/write_description.R
	Rscript dev/install_playbase.R 'github'

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
	docker build $(ARG) --no-cache -f dev/Dockerfile.os \
	    -t playbase-os . 2>&1 | tee docker-os.log

docker.rbase: 
	docker build $(ARG) --no-cache -f dev/Dockerfile.rbase \
	    -t playbase-rbase . 2>&1 | tee docker-rbase.log

docker.pkg: 
	docker build $(ARG) -f dev/Dockerfile \
	    -t playbase-pkg . 2>&1 | tee docker.log

## we need to squash the layer to minimize the size but also to hide
## the use of any temporary tokens. Install docker-squash from pipx.
docker.squash: 
	@if [ -z `command -v pipx &> /dev/null` ]; then \
	    echo ERROR: please install docker-squash; \
	    exit 1; \
	fi
	docker-squash playbase-pkg -t bigomics/playbase:latest

docker: docker.os docker.rbase docker.pkg docker.squash
	@echo "% building playbase docker..."
