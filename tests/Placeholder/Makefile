ROOT_DIR =	../..
INCS	=	-I$(ROOT_DIR)/src
DEFS    =	-g -O3 -Wno-deprecated -Wall -std=c++17 -fopenmp -DMPIPARALLEL #
LIBS	=	-lm
COMP	=	mpicxx
DEPS    =       
OBJ     =	main.cc


all: run2.exe 

%.o: %.cc $(DEPS)
	$(COMP) -c -o $@ $< $(DEFS) $(LIBS)

run2.exe: $(OBJ)
	$(COMP) -o run2.exe $(OBJ) $(INCS) $(DEFS) $(LIBS)

clean:
	rm -f *.o run2.exe
	

#rm -r data