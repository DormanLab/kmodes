# 'Makefile'
MARKDOWN = pandoc --from gfm --to html --standalone --metadata pagetitle="k-modes"
all: $(patsubst %.md,%.html,$(wildcard *.md)) Makefile

clean:
	rm -f $(patsubst %.md,%.html,$(wildcard *.md))
	rm -f *.bak *~

%.html: %.md
	$(MARKDOWN) $< --output $@
