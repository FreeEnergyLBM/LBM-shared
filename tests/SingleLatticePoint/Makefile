
INCS	=	

DEFS    =	   -O2 -Wno-deprecated -Wno-unused-result -Wall -std=c++17  -I.  

LIBS	=	-lm

LIBPATH =

COMP	=	 mpicxx #mpic++ 

DEPS    =       ../../src/%.hh ../../src/Forces/%.hh ../../src/GradientStencils%.hh ../../src/LBModels/%.hh ../../src/BoundaryModels/%.hh 

OBJ     =       main.cc

all: wetP 
#all: wetS 

%.o: %.cc $(DEPS)
	$(COMP) -c -o $@ $< $(DEFS) $(LIBS)

wetP: $(OBJ)
	$(COMP) -o wetP $(OBJ)  $(INCS) $(DEFS) $(LIBPATH) $(LIBS)

#wetS: $(OBJ)
#	$(COMP) -o wetS $(OBJ)  $(INCS) $(DEFS) $(LIBPATH) $(LIBS)

clean:
	rm -f *.o  wetP wetS wetS.exe *~ *.bak	

