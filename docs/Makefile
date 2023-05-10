all: html latex

html:
	(cat Doxyfile ; echo "GENERATE_HTML=YES") | doxygen -

latex:
	(cat Doxyfile ; echo "GENERATE_LATEX=YES") | doxygen -
	$(MAKE) -C latex
	cp latex/refman.pdf doxygen.pdf

clean:
	rm -rf html latex doxygen.pdf
