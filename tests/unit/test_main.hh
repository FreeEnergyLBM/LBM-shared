#include <vector>

#include "Mpi.hh"
#include "gtest-mpi-listener.hpp"
#include "gtest/gtest.h"

template <typename T>
::testing::AssertionResult ArraysMatch(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size()) {
        return ::testing::AssertionFailure() << "array 1 has size " << a.size() << " != expected size " << b.size();
    }
    if (a != b) {
        return ::testing::AssertionFailure() << ::testing::PrintToString(a) << " != " << ::testing::PrintToString(b);
    }
    return ::testing::AssertionSuccess();
}

template <typename T>
::testing::AssertionResult ArraysNear(const std::vector<T>& a, const std::vector<T>& b, float delta = 1e-6) {
    if (a.size() != b.size()) {
        return ::testing::AssertionFailure() << "array 1 has size " << a.size() << " != expected size " << b.size();
    }
    for (size_t i(0); i < a.size(); ++i) {
        if (std::abs(a[i] - b[i]) < delta) continue;
        return ::testing::AssertionFailure() << ::testing::PrintToString(a) << " != " << ::testing::PrintToString(b);
    }
    return ::testing::AssertionSuccess();
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

#ifdef MPIPARALLEL
    MPI_Init(&argc, &argv);
    mpi.init(MPI_COMM_WORLD);

    // Add an MPI listener (https://github.com/LLNL/gtest-mpi-listener)
    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    ::testing::TestEventListener* l = listeners.Release(listeners.default_result_printer());
    listeners.Append(new GTestMPIListener::MPIWrapperPrinter(l, MPI_COMM_WORLD));
#endif

    return RUN_ALL_TESTS();
}
