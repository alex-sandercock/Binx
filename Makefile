BINX_OUT ?= binx_gwas.tsv
GWASPOLY_OUT ?= gwaspoly_results.tsv
TRAIT ?= Trait
OUT ?= benchmark/merged.tsv

.PHONY: benchmark-gwaspoly
benchmark-gwaspoly:
	@echo "Comparing binx ($(BINX_OUT)) vs GWASpoly ($(GWASPOLY_OUT))"
	@Rscript scripts/compare_binx_gwaspoly.R --binx=$(BINX_OUT) --gwaspoly=$(GWASPOLY_OUT) --trait=$(TRAIT) --out=$(OUT)
