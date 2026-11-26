PYTHON ?= python3

.PHONY: all results figures article pdf test clean

all: pdf

results:
	$(PYTHON) code/calculs.py code data results

figures: results
	$(PYTHON) code/figures.py results paper/figures

article: figures
	$(PYTHON) code/remplir.py results/resultats.json paper/parts code/velib.py paper/article.tex

pdf: article
	cd paper && pdflatex -interaction=nonstopmode article.tex && pdflatex -interaction=nonstopmode article.tex

test:
	$(PYTHON) -m pytest -q

clean:
	rm -f results/series.npz paper/*.aux paper/*.log paper/*.out
