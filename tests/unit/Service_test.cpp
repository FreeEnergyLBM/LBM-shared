// YPROCS 2
#include "Service.hh"

#include "Lattice.hh"
#include "Mpi.hh"
#include "test_main.hh"

using Lattice = LatticeProperties<ParallelY<>, 2, 4, 3>;

TEST(Service, computeXYZ) {
    Lattice();

    std::array<int, 3> xyz;
    xyz = computeXYZ<Lattice>(0);
    EXPECT_EQ(xyz[0], 0);
    EXPECT_EQ(xyz[1], (mpi.rank == 0) ? 3 : 1);
    EXPECT_EQ(xyz[2], 0);

    xyz = computeXYZ<Lattice>(4);
    EXPECT_EQ(xyz[0], 0);
    EXPECT_EQ(xyz[1], (mpi.rank == 0) ? 0 : 2);
    EXPECT_EQ(xyz[2], 1);

    xyz = computeXYZ<Lattice>(12);
    EXPECT_EQ(xyz[0], 1);
    EXPECT_EQ(xyz[1], (mpi.rank == 0) ? 3 : 1);
    EXPECT_EQ(xyz[2], 0);
}

TEST(Service, computeK) {
    Lattice();

    EXPECT_EQ(computeK<Lattice>(0, 0, 0), 3);
    EXPECT_EQ(computeK<Lattice>(0, 1, 1), 7);
    EXPECT_EQ(computeK<Lattice>(1, 1, 2), 20);
}

TEST(Service, RangeIterator) {
    Lattice();

    std::vector<int> xs, ys, zs, ks;
    for (auto [x, y, z, k] : RangeXYZK<Lattice>()) {
        xs.push_back(x);
        ys.push_back(y);
        zs.push_back(z);
        ks.push_back(k);
    }
    auto ys_test = (mpi.rank == 0) ? std::vector<int>{0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1}
                                   : std::vector<int>{2, 2, 2, 3, 3, 3, 2, 2, 2, 3, 3, 3};
    EXPECT_TRUE(ArraysMatch(xs, {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}));
    EXPECT_TRUE(ArraysMatch(ys, ys_test));
    EXPECT_TRUE(ArraysMatch(zs, {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2}));
    EXPECT_TRUE(ArraysMatch(ks, {3, 4, 5, 6, 7, 8, 15, 16, 17, 18, 19, 20}));
}

template <int i>
struct TestFunction {
    static bool called;
    static void call() { called = true; }
};

template <int i>
bool TestFunction<i>::called = false;

TEST(Service, constexpr_for) {
    constexpr_for<3>([]<int i>() { TestFunction<i>::call(); });
    EXPECT_TRUE(TestFunction<0>::called);
    EXPECT_TRUE(TestFunction<1>::called);
    EXPECT_TRUE(TestFunction<2>::called);
    EXPECT_FALSE(TestFunction<3>::called);
}
