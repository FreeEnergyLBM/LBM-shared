ROOT_DIR = ../..
INC_DIR = $(ROOT_DIR)/src
INC = $(addprefix -I, $(INC_DIR))

TESTS = ParallelX_test ParallelY_test Collide_test BounceBack_test CubicWetting_test FlowField_test
RUN_TESTS = $(addprefix run_, $(TESTS))

CXX = mpicxx
CXXFLAGS = -std=c++17 -Wall -DMPIPARALLEL -fopenmp

# Read the number of processors from the top of the file. Defaults to 1
define GETPROCS
$(eval XPROCS=$(shell head -n 1 $(1) | sed -n 's/.*XPROCS \([0-9]\+\).*/\1/p'))
$(eval YPROCS=$(shell head -n 1 $(1) | sed -n 's/.*YPROCS \([0-9]\+\).*/\1/p'))
$(eval ZPROCS=$(shell head -n 1 $(1) | sed -n 's/.*ZPROCS \([0-9]\+\).*/\1/p'))
$(if $(XPROCS), , $(eval XPROCS=1))
$(if $(YPROCS), , $(eval YPROCS=1))
$(if $(ZPROCS), , $(eval ZPROCS=1))
$(eval NPROCS = $(shell echo $$(($(XPROCS)*$(YPROCS)*$(ZPROCS)))))
endef

.PHONY: all clean cleantests gtest $(RUN_TESTS)

# Use 'failed' file to determine success so that all tests are performed even if one fails
all: $(TESTS) $(RUN_TESTS) cleantests
	@[ -f failed ] && rm failed && exit 1 || exit 0

clean:
	rm -f $(TESTS) gtest.a gtest_main.a *.o

cleantests:
	rm -f $(TESTS) *_test.o

$(RUN_TESTS): run_%: %
	$(call GETPROCS, $*.cpp)
	@echo "Running $* on $(NPROCS) process(es)"
	@mpirun -np $(NPROCS) ./$* || touch failed

# Googletest rules and parameters, see:
# https://github.com/google/googletest/blob/release-1.8.1/googletest/make/Makefile
GTEST_DIR = googletest/googletest
CPPFLAGS += -isystem $(GTEST_DIR)/include
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
# Googletest with MPI
# https://github.com/LLNL/gtest-mpi-listener
GTEST_MPI_DIR = gtest-mpi-listener
INC += -I$(GTEST_MPI_DIR)/include

gtest:
	git submodule update --init $(dir $(GTEST_DIR))
	git submodule update --init $(GTEST_MPI_DIR)

gtest-all.o: $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o: $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest_main.cc

gtest.a: gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a: gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Tests
$(TESTS): %_test: %_test.o gtest.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -lpthread -o $@

%_test.o: %_test.cpp $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INC) -c $<
