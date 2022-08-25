SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

RSTW = $(shell find . -name "*rstw")
RST = $(RSTW:.rstw=.rst)

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Avoid catching these files in Sphinx rules
%.rstw:

%.rst: %.rstw
	pweave -F $(shell dirname $<) $< -o $@

%: Makefile $(RST)
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
