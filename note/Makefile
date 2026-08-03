# Build the note.  Figures come from ../docs/figures, which are produced by
#   python3 python/make_xreco_report.py --png-dir docs/figures
# so regenerate those first if the numbers have changed.
all: fasercal_muondis.pdf

fasercal_muondis.pdf: fasercal_muondis.tex
	pdflatex -interaction=nonstopmode $< >/dev/null
	pdflatex -interaction=nonstopmode $< >/dev/null

clean:
	rm -f *.aux *.log *.out *.toc
