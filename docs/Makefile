.PHONY: all clean doxygen figures html

FIGURES = rare.png isochromat.png epg/regular.png epg/discrete.png epg/discrete_3d.png

all: html

clean:
	rm -rf _build $(FIGURES)

doxygen:
	mkdir -p _build
	doxygen

figures: $(FIGURES)

html: doxygen $(FIGURES)
	mkdir -p _build
	python3 -m sphinx -T -E -b html -d _build/doctrees . _build/html

epg/%.png: ../examples/epg/%.py
	python3 $< $@

%.png: ../examples/%.py
	python3 $< $@
