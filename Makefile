.PHONY: help install trf1 trf1-no-chrM clean

CONDA_ENV ?= Bioinfomatics
PY ?= conda run -n $(CONDA_ENV) python
PEAK := data/GSM638201_TRF1-p0.05.bed.gz
DEG := data/siTRF1_DEG_list_gene_name.txt
ANNOTATION := data/gene_annotation.tsv
OUTDIR ?= results

help:
	@echo "Targets:"
	@echo "  install        pip install -e . into the active environment"
	@echo "  trf1           Run the full TRF1 promoter overlap pipeline -> $(OUTDIR)/"
	@echo "  trf1-no-chrM   Same as trf1 but excludes mitochondrial peaks -> results_no_chrM/"
	@echo "  clean          Remove $(OUTDIR)/ and results_no_chrM/"

install:
	$(PY) -m pip install -e .

trf1:
	$(PY) -m trf1_overlap \
	    --peak $(PEAK) \
	    --deg $(DEG) \
	    --annotation $(ANNOTATION) \
	    --outdir $(OUTDIR)

trf1-no-chrM:
	$(PY) -m trf1_overlap \
	    --peak $(PEAK) \
	    --deg $(DEG) \
	    --annotation $(ANNOTATION) \
	    --outdir results_no_chrM \
	    --exclude-chroms chrM

clean:
	rm -rf $(OUTDIR) results_no_chrM
