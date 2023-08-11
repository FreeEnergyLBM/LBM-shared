ROOT_DIR = ../../../..
INC_DIR = $(ROOT_DIR)/src
INC = $(addprefix -I, $(INC_DIR))

CXX = mpicxx
CXXFLAGS = -O3 -std=c++17 -Wall -DMPIPARALLEL -fopenmp -DMPIPARALLEL

all: run.exe

run.exe: main.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INC) $^ -lpthread -o $@

clean:
	rm -f run.exe
