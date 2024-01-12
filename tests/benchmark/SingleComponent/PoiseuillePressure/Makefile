ROOT_DIR = ../../../..
INC_DIR = $(ROOT_DIR)/src
INC = $(addprefix -I, $(INC_DIR))

CXX = mpicxx
CXXFLAGS = -O3 -std=c++17 -Wall -DMPIPARALLEL -fopenmp -DMPIPARALLEL

ifdef strip_templates
STRIP_TEMPLATES = -fdiagnostics-color=always 2>&1 | sed 's/ \[with [^]]*\]//g'
endif

run.exe: main.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INC) $^ -lpthread -o $@ $(STRIP_TEMPLATES)

clean:
	rm -f run.exe
	rm -r data
