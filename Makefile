.PHONY: all clean

all: README.md vignettes/examples.md vignettes/examples.html

clean:
	rm -rf README_{files,cache} README.{md,html} vignettes/examples_* vignettes/examples.{md,html}

README.md: README.Rmd
	R -e 'rmarkdown::render("README.Rmd", "all")'

vignettes/examples.md: vignettes/examples.Rmd
	R -e 'devtools::build_rmd("vignettes/examples.Rmd")'

vignettes/examples.html: vignettes/examples.Rmd vignettes/custom.css
	R -e 'devtools::build_rmd("vignettes/examples.Rmd")'
