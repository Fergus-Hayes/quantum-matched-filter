DOCS = quantum.pdf

PACK_FILES = quantum.tex $(FIGS)

FIGS = 

TABLES = 

BIBLIO = masterbib.bib

DEPS = $(FIGS) $(TABLES) $(BIBLIO)

BYPROD = *.aux *.log *.out *.dvi *.blg *.bbl __tmp \
         quantum.bbl quantum.blg quantum.ps \
         quantumNotes.bib tag.tex quantum.tex quantum.pdf __tag.tex

all : git-tag $(DOCS)

.PHONY : git-tag

default : all

git-tag :
	./git-tag.sh

quantum.tex: quantum.tex quantum.bbl $(TABLES)
	latexpand --output paper.tex --expand-bbl quantum.bbl quantum.tex

quantum.pdf quantum.bbl: quantum.tex $(DEPS)
	pdflatex quantum && bibtex quantum && pdflatex quantum && pdflatex quantum && pdflatex quantum

clean :
	rm -f $(BYPROD) $(DOCS)
