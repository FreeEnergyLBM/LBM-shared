ROOT_DIR = ../..
INC_DIR = $(ROOT_DIR)/src
INC = $(addprefix -I, $(INC_DIR))

SHELL = /bin/bash -O globstar
CXX = mpicxx
CXXFLAGS = -std=c++17 -Wall -DMPIPARALLEL -fopenmp

BENCHMARKS = SingleComponent/PoiseuilleForce SingleComponent/PoiseuillePressure Binary/LaplacePressure Binary/LayeredPoiseuille

.PHONY: all clean cleantests gtest $(BENCHMARKS)

# Use 'failed' file to determine success so that all tests are performed even if one fails
all: $(BENCHMARKS) clean
	@[ -f failed ] && rm failed && exit 1 || exit 0

clean:
	rm -f **/run.exe failed

$(BENCHMARKS): %: %/run.exe
	@echo "Running $*"
	@mpirun -np 1 $< && python3 $*/Analysis.py && (echo -e "\033[0;32mSUCCESS\033[0m") || (echo -e "\033[0;31mFAILED\033[0m" && touch failed)

BENCH_RUN = $(addsuffix /run.exe, $(BENCHMARKS))
$(BENCH_RUN): %/run.exe: %/main.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INC) $^ -lpthread -o $@
