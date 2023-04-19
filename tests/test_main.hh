#include "gtest/gtest.h"
#include "gtest-mpi-listener.hpp"


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);

  #ifdef MPIPARALLEL
  MPI_Init(&argc, &argv);

  // Add an MPI listener (https://github.com/LLNL/gtest-mpi-listener)
  ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  ::testing::TestEventListener *l = listeners.Release(listeners.default_result_printer());
  listeners.Append(new GTestMPIListener::MPIWrapperPrinter(l, MPI_COMM_WORLD));
  #endif

  return RUN_ALL_TESTS();
}
