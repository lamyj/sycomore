.PHONY: clean doxygen html

FIGURES = epg/regular.png epg/discrete.png epg/discrete_3d.png

all: html

clean:
	rm -rf _build $(FIGURES)

doxygen:
	mkdir -p _build
	doxygen

html: doxygen $(FIGURES)
	mkdir -p _build
	python3 -m sphinx -T -E -b html -d _build/doctrees . _build/html

%.png: %.py
	python3 $< $@
