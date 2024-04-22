
.PHONY: clean all

all: main.pdf

main.pdf: main.tex
	for i in 1 2; \
	do \
		pdflatex -output-directory=pdflatex -interaction=nonstopmode main.tex &> pdflatex/stderr.txt; \
	done; mv pdflatex/main.pdf ./

clean:
	command rm -f main.pdf
	command rm -rf pdflatex/*
