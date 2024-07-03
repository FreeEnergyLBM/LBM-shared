ROOT_DIR =	../..
INCS	=	-I$(ROOT_DIR)/src
DEFS    =	-g -O3 -Wno-deprecated -Wall -fmax-errors=1 -std=c++17 -fopenmp -DMPIPARALLEL
LIBS	=	-lm
COMP	=	mpicxx
DEPS    =       
OBJ     =	main.cc

ifdef strip_templates
STRIP_TEMPLATES = -fdiagnostics-color=always 2>&1 | sed 's/ \[with [^]]*\]//g'
endif


all: run.exe 

%.o: %.cc $(DEPS)
	$(COMP) -c -o $@ $< $(DEFS) $(LIBS)

run.exe: $(OBJ)
	$(COMP) -o run.exe $(OBJ) $(INCS) $(DEFS) $(LIBS) $(STRIP_TEMPLATES)

clean:
	rm -f *.o run.exe
	rm -r data

