# makefile to generate documentation
# Author: Cheng Li

SCRIPT_DIR = scripts/
EXAMPLE_DIR = examples/
PGEN_DIR = ../src/pgen/
DOXYGEN_DIR = doxygen/
COUNT = 1

EXAMPLES := 2d.straka
EXAMPLE_DOCS := $(addprefix $(EXAMPLE_DIR),$(EXAMPLES:%=%.dox))
PLAIN_DOCS := $(patsubst 0d.%,%,$(EXAMPLES))
PLAIN_DOCS := $(patsubst 1d.%,%,$(PLAIN_DOCS))
PLAIN_DOCS := $(patsubst 2d.%,%,$(PLAIN_DOCS))
PLAIN_DOCS := $(patsubst 3d.%,%,$(PLAIN_DOCS))
PLAIN_DOCS := $(addsuffix _plain.txt,$(addprefix $(EXAMPLE_DIR),$(PLAIN_DOCS)))

all : $(EXAMPLE_DOCS) $(PLAIN_DOCS) doxygen

$(EXAMPLE_DIR)%.dox: ../$(EXAMPLE_DIR)%
	mkdir -p $(EXAMPLE_DIR)
	perl $(SCRIPT_DIR)make_step.pl $(subst ../$(EXAMPLE_DIR),,$<) ../ $(COUNT) > $@
	@ COUNT=$(shell expr $(COUNT) + 1)
	@ mkdir -p $(DOXYGEN_DIR)img
	@ cp $(<)/*.png $(DOXYGEN_DIR)img

$(EXAMPLE_DIR)%_plain.txt: $(PGEN_DIR)%.cpp
	mkdir -p $(EXAMPLE_DIR)
	perl $(SCRIPT_DIR)program2plain $< > $@

doxygen: $(DOXYGEN_DIR)configure
	cd doxygen && doxygen configure
	@ cp -r $(DOXYGEN_DIR)img $(DOXYGEN_DIR)html

.PHONY : all doxygen

clean:
	rm -f $(EXAMPLE_DIR)*.dox
	rm -f $(EXAMPLE_DIR)*_plain.txt
	rm -f $(DOXYGEN_DIR)bib*.aux
	rm -f $(DOXYGEN_DIR)citelist.doc*
	rm -rf $(DOXYGEN_DIR)html
	rm -rf $(DOXYGEN_DIR)latex
	rm -rf $(DOXYGEN_DIR)img
