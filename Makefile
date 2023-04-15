build: doc
	R -e "devtools::build()"

doc:
	R -e "devtools::document()"
	R -e "devtools::build_vignettes()"

check: 
	R -e "devtools::check()"

install:
	R CMD INSTALL .
#	R CMD INSTALL ../playbase_0.1.0.tar.gz

