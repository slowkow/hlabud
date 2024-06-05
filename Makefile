.PHONY: all clean install check test

# ROOT_DIR:=$(shell dirname "$(realpath $(firstword $(MAKEFILE_LIST)))")

rmds := $(wildcard vignettes/articles/*.Rmd)
articles = $(rmds:.Rmd=.html)
css = vignettes/articles/custom.css
mans := $(wildcard man/*.Rd)

all: README.md docs/index.html

clean:
	rm -rf docs/ README_{files,cache} README.{md,html} index.{md,html} vignettes/articles/examples_{files,cache} vignettes/articles/*.html

install:
	R -e 'devtools::document(); devtools::install()'

check:
	R -e 'devtools::document(); rcmdcheck::rcmdcheck(args="--no-tests")'

cran:
	R -e 'devtools::document(); rcmdcheck::rcmdcheck(args="--as-cran")'

test:
	R -e 'devtools::test()'

README.md: README.Rmd $(articles)
	# R -e 'devtools::install_deps(".", TRUE)'
	R -e 'rmarkdown::render("README.Rmd", "all")'

index.md: index.Rmd $(articles)
	R -e 'rmarkdown::render("index.Rmd", "md_document")'

vignettes/articles/%.html: vignettes/articles/%.Rmd $(css)
	R -e 'devtools::build_rmd("$<")'

# vignettes/articles/examples.html: vignettes/articles/examples.Rmd vignettes/articles/custom.css
# 	R -e 'devtools::build_rmd("vignettes/articles/examples.Rmd")'

# vignettes/articles/visualize-hla-structure.html: vignettes/articles/visualize-hla-structure.Rmd
# 	R -e 'devtools::build_rmd("vignettes/articles/visualize-hla-structure.Rmd")'

# vignettes/articles/numbering.html: vignettes/articles/numbering.Rmd
# 	R -e 'devtools::build_rmd("vignettes/articles/numbering.Rmd")'

docs/index.html: index.md $(mans)
	rm -rf docs/
	R -e 'pkgdown::init_site(); pkgdown::build_articles(); pkgdown::build_site()'
	rm -f docs/paper.*

