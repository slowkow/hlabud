.PHONY: all clean install

all: README.md vignettes/examples.md vignettes/examples.html docs/articles/examples.html

clean:
	rm -rf README_{files,cache} README.{md,html} vignettes/examples_{files,cache} vignettes/examples.{md,html} docs/

install:
	R -e 'devtools::document(); devtools::install()'

README.md: README.Rmd
	R -e 'devtools::install_deps(".", TRUE)'
	R -e 'rmarkdown::render("README.Rmd", "all")'

vignettes/examples.md: vignettes/examples.Rmd
	# R -e 'devtools::build_rmd("vignettes/examples.Rmd", ")'
	R -e 'rmarkdown::render("vignettes/examples.Rmd", output_format = "md_document")'

vignettes/examples.html: vignettes/examples.Rmd vignettes/custom.css
	R -e 'devtools::build_rmd("vignettes/examples.Rmd")'

docs/articles/examples.html: vignettes/examples.Rmd index.md man/*.Rd
	R -e 'pkgdown::build_site()'

