ROOT_DIR =       ../..

INCS	=	-I$(ROOT_DIR)/src

DEFS    =	-g -O3 -Wno-deprecated -Wall -std=c++17 -fopenmp -DMPIPARALLEL

LIBS	=	-lm

COMP	=	mpicxx

DEPS    =       

OBJ     =	main.cc


all: run.exe 

%.o: %.cc $(DEPS)
	$(COMP) -c -o $@ $< $(DEFS) $(LIBS)

run.exe: $(OBJ)
	$(COMP) -o run.exe $(OBJ) $(INCS) $(DEFS) $(LIBS)

clean:
	rm -f *.o run.exe
	rm -r data

