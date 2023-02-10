
INCS	=	

DEFS    =	   -O2 -Wno-deprecated -Wno-unused-result -Wall -std=c++17  -I.  

LIBS	=	-lm

LIBPATH =

COMP	=	 mpicxx #mpic++ 

DEPS    =       ../../src/Algorithm.hh ../../src/Collide.hh ../../src/Data.hh ../../src/Forces.hh ../../src/Geometry.hh ../../src/Global.hh ../../src/LBModels.hh ../../src/Parameter.hh ../../src/Service.hh ../../src/Stencil.hh # RT_Timer.hh Timer.hh Bool.hh String.hh 

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

