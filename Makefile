.PHONY: all clean install check test

# ROOT_DIR:=$(shell dirname "$(realpath $(firstword $(MAKEFILE_LIST)))")

all: README.md vignettes/articles/examples.md vignettes/articles/examples.html docs/articles/examples.html vignettes/articles/visualize-hla-structure.html

clean:
	rm -rf README_{files,cache} README.{md,html} vignettes/articles/examples_{files,cache} vignettes/articles/examples.{md,html} docs/

install:
	R -e 'devtools::document(); devtools::install()'

check:
	R -e 'devtools::document(); rcmdcheck::rcmdcheck(args="--no-tests")'

test:
	R -e 'devtools::test()'

README.md: README.Rmd
	R -e 'devtools::install_deps(".", TRUE)'
	R -e 'rmarkdown::render("README.Rmd", "all")'

vignettes/articles/examples.md: vignettes/articles/examples.Rmd
	# R -e 'devtools::build_rmd("vignettes/articles/examples.Rmd", ")'
	R -e 'rmarkdown::render("vignettes/articles/examples.Rmd", output_format = "md_document")'

vignettes/articles/examples.html: vignettes/articles/examples.Rmd vignettes/articles/custom.css
	R -e 'devtools::build_rmd("vignettes/articles/examples.Rmd")'

docs/articles/examples.html: vignettes/articles/examples.Rmd index.md man/*.Rd
	R -e 'pkgdown::build_site()'
	rm -f docs/paper.*

vignettes/articles/visualize-hla-structure.html: vignettes/articles/visualize-hla-structure.Rmd
	R -e 'devtools::build_rmd("vignettes/articles/visualize-hla-structure.Rmd")'
