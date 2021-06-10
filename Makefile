# 'Makefile'
MARKDOWN = pandoc --from gfm --to html --standalone --metadata pagetitle="k-modes"
all: $(patsubst %.md,%.html,$(wildcard *.md)) Makefile

clean:
	rm -f $(patsubst %.md,%.html,$(wildcard *.md))
	rm -f *.bak *~
	rm -f src/kmodes.so src/*.o src/cmake_install.cmake src/CMakeCache.txt src/run_kmodes src/Makefile
	rm -fr src/CMakeFiles

%.html: %.md
	$(MARKDOWN) $< --output $@
