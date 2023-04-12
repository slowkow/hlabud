.PHONY: all clean

all: README.md

clean:
	rm -rf README_files README_cache README.md README.html

README.md: README.Rmd
	R -e 'rmarkdown::render("README.Rmd", "all")'

