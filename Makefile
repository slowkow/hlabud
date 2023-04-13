.PHONY: all clean

all: README.md vignettes/examples.md vignettes/examples.html docs/articles/examples.html

clean:
	rm -rf README_{files,cache} README.{md,html} vignettes/examples_* vignettes/examples.{md,html} docs/

README.md: README.Rmd
	R -e 'rmarkdown::render("README.Rmd", "all")'

vignettes/examples.md: vignettes/examples.Rmd
	R -e 'devtools::build_rmd("vignettes/examples.Rmd")'

vignettes/examples.html: vignettes/examples.Rmd vignettes/custom.css
	R -e 'devtools::build_rmd("vignettes/examples.Rmd")'

docs/articles/examples.html: vignettes/examples.Rmd index.md man/*.Rd
	R -e 'pkgdown::build_site()'

