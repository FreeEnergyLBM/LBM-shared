// XPROCS 4
#include "Lattice.hh"
#include "Parallel.hh"
#include "test_main.hh"

constexpr int XSIZE = 128;  // global size along X, parallelised
constexpr int ysize = 128;  // global==local since no parallelisation along Y
constexpr int zsize = 128;  // global==local since no parallelisation along Z
constexpr int width = 2;
using Parallel_Pattern = ParallelX<width>;
using Lattice = LatticeProperties<Parallel_Pattern, XSIZE, ysize, zsize>;

// In the tests below:
// - I introduce the xExt variable ('x extended')
//   since parallelisation is done along X and the local size along X includes 'width' from both sides.
// - 'ysize' and 'zsize' remain the same full sizes since parallelisation isn't done along Y and Z.

// In the 'TestCommunicateParameter' TEST:
// - I construct the 2 vectors:
//   (1) for communication ('density') and
//   (2) for testing against (1) ('densityArray' containing the expected values).
// - The communication vectors are constructed such that its elements have consecutive numbers beginning from 1
//   (ordered firstly along Z, secondly along Y, thirdly along X)

TEST(ParallelX_Test, TestCommunicateParameter) {
    Lattice::ResetParallelTracking();

    Lattice lattice;
    std::vector<double> &density = Density<>::getInstance<Lattice>().mv_Parameter;
    int xExt = Lattice::LXdiv;
    ASSERT_EQ(ysize, Lattice::LYdiv);
    ASSERT_EQ(zsize, Lattice::LZdiv);
    ASSERT_EQ(density.size(), xExt * ysize * zsize);
    std::fill(density.begin(), density.end(), 0.0);
    for (int x = width; x < xExt - width; ++x)
        for (int y = 0; y < ysize; ++y)
            for (int z = 0; z < zsize; ++z) {
                int i = z + y * zsize + x * ysize * zsize;
                // rank=0:  z + y*zsize + (x-width)*ysize*zsize + 1
                // rank>0:  z + y*zsize + (x-width)*ysize*zsize + (mpi.rank*(xExt-2*width))*ysize*zsize + 1
                density[i] = i + (mpi.rank * (xExt - 2 * width) - width) * ysize * zsize + 1;
                //          z + y*zsize + x*ysize*zsize + (mpi.rank*(xExt-2*width)-width)*ysize*zsize + 1
            }
    std::vector<double> densityArray(density.begin(), density.end());
    for (int x = 0; x < width; ++x)
        for (int y = 0; y < ysize; ++y)
            for (int z = 0; z < zsize; ++z) {
                int i = z + y * zsize + x * ysize * zsize;
                densityArray[i] = density[i + width * ysize * zsize];
                densityArray[i] += (mpi.rank > 0) ? -width * ysize * zsize : (XSIZE - width) * ysize * zsize;
            }
    for (int x = xExt - width; x < xExt; ++x)
        for (int y = 0; y < ysize; ++y)
            for (int z = 0; z < zsize; ++z) {
                int i = z + y * zsize + x * ysize * zsize;
                densityArray[i] = density[i - width * ysize * zsize];
                densityArray[i] += (mpi.rank < mpi.size - 1) ? width * ysize * zsize : (width - XSIZE) * ysize * zsize;
            }

    lattice.communicate(Density<>::getInstance<Lattice>());
    EXPECT_TRUE(ArraysMatch(density, densityArray));
}

// In the 'TestCommunicateDistribution*' TESTs:
// - I construct the 2 sets of vectors:
//   (1) for communication ('distrNodeLeft' and 'distrNodeRight') and
//   (2) for testing against (1) ('nodeLeftArray' and 'nodeRightArray' containing the expected values).
// - The communication vectors are constructed such that its elements have consecutive numbers beginning from 0
//   (ordered firstly along velocity directions, and then along Z, along Y, and along X)

//
// D2Q5
//
TEST(ParallelX_Test, TestCommunicateDistributionD2Q5) {
    Lattice::ResetParallelTracking();

    Lattice lattice;
    using DQ = D2Q5;
    constexpr int Q = DQ::Q;
    constexpr int neighbors = Parallel_Pattern::mNumDirections * 2;
    DataOldNew<Lattice, DQ> data;
    int xExt = Lattice::LXdiv;  // length of local buffer including external halo region along X
    int nodeLeft = (width - 1) * ysize * zsize * Q;
    int nodeRight = (xExt - width) * ysize * zsize * Q;

    auto distr = data.getDistributionObject();
    ASSERT_EQ(ysize, Lattice::LYdiv);
    ASSERT_EQ(zsize, Lattice::LZdiv);
    ASSERT_EQ(distr.mv_Distribution.size(), xExt * ysize * zsize * Q);  // #elements in the distribution
    std::fill(distr.mv_Distribution.begin(), distr.mv_Distribution.end(), 0.0);

    // 1. DISTRIBUITION INITIALISATION FOR COMMUNICATION
    for (int i = 0; i < neighbors; ++i)
        for (int y = 0; y < ysize; ++y)
            for (int z = 0; z < zsize; ++z)
                for (int j = 0; j < Q; ++j)  // #directions in the stencil
                {
                    int k = (i == 0)
                                ? nodeLeft + (z + y * zsize) * Q + j
                                : nodeRight + (z + y * zsize) * Q + j;  // element in the distribution (with directions)
                    double value = ((neighbors * mpi.rank + i) * ysize * zsize + (z + y * zsize)) * Q +
                                   j;  // the value in the element of the distribution
                    distr.getDistribution(k) = value;
                }

    // 2. SETTING THE EXPECTED VECTOR WITH EXPECTED VALUES AFTER COMMUNICATION
    std::vector<double> nodeLeftArray(Q * ysize * zsize, 0.0);
    std::vector<double> nodeRightArray(Q * ysize * zsize, 0.0);
    for (int y = 0; y < ysize; ++y)
        for (int z = 0; z < zsize; ++z) {
            int k = (z + y * zsize) * Q;
            nodeLeftArray[k + 1] = distr.getDistribution(nodeLeft + k + 1);
            nodeLeftArray[k + 1] += (mpi.rank > 0) ? -ysize * zsize * Q : (2 * mpi.size - 1) * ysize * zsize * Q;
            nodeRightArray[k + 2] = distr.getDistribution(nodeRight + k + 2);
            nodeRightArray[k + 2] +=
                (mpi.rank < mpi.size - 1) ? ysize * zsize * Q : (1 - 2 * mpi.size) * ysize * zsize * Q;
        }

    // 3. DISTRIBUITION COMMUNICATION
    lattice.communicateDistribution(distr);

    // 4. RETRIEVING DATA FROM DISTRIBUTION AFTER COMMUNICATION
    std::vector<double> distrNodeLeft(distr.getDistributionPointer(width * ysize * zsize),
                                      distr.getDistributionPointer((width + 1) * ysize * zsize));
    std::vector<double> distrNodeRight(distr.getDistributionPointer((xExt - width - 1) * ysize * zsize),
                                       distr.getDistributionPointer((xExt - width) * ysize * zsize));

    // 5. COMPARISON WITH THE EXPECTED VALUES
    EXPECT_TRUE(ArraysMatch(distrNodeLeft, nodeLeftArray));
    EXPECT_TRUE(ArraysMatch(distrNodeRight, nodeRightArray));
}

//
// D2Q9
//
TEST(ParallelX_Test, TestCommunicateDistributionD2Q9) {
    Lattice::ResetParallelTracking();

    Lattice lattice;
    using DQ = D2Q9;
    constexpr int Q = DQ::Q;
    constexpr int neighbors = Parallel_Pattern::mNumDirections * 2;
    DataOldNew<Lattice, DQ> data;
    int xExt = Lattice::LXdiv;  // length of local buffer including external halo region along X
    int nodeLeft = (width - 1) * ysize * zsize * Q;
    int nodeRight = (xExt - width) * ysize * zsize * Q;

    auto distr = data.getDistributionObject();
    ASSERT_EQ(ysize, Lattice::LYdiv);
    ASSERT_EQ(zsize, Lattice::LZdiv);
    ASSERT_EQ(distr.mv_Distribution.size(), xExt * ysize * zsize * Q);  // #elements in the distribution
    std::fill(distr.mv_Distribution.begin(), distr.mv_Distribution.end(), 0.0);

    // 1. DISTRIBUITION INITIALISATION FOR COMMUNICATION
    for (int i = 0; i < neighbors; ++i)
        for (int y = 0; y < ysize; ++y)  // #elements along Y axis
            for (int z = 0; z < zsize; ++z)
                for (int j = 0; j < Q; ++j)  // #directions in the stencil
                {
                    int k = (i == 0)
                                ? nodeLeft + (z + y * zsize) * Q + j
                                : nodeRight + (z + y * zsize) * Q + j;  // element in the distribution (with directions)
                    double value = ((neighbors * mpi.rank + i) * ysize * zsize + (z + y * zsize)) * Q +
                                   j;  // the value in the element of the distribution
                    distr.getDistribution(k) = value;
                }

    // 2. SETTING THE EXPECTED VECTOR WITH EXPECTED VALUES AFTER COMMUNICATION
    std::vector<double> nodeLeftArray(Q * ysize * zsize, 0.0);
    std::vector<double> nodeRightArray(Q * ysize * zsize, 0.0);
    for (int y = 0; y < ysize; ++y)
        for (int z = 0; z < zsize; ++z) {
            int k = (z + y * zsize) * Q;
            nodeLeftArray[k + 1] = distr.getDistribution(nodeLeft + k + 1);
            nodeLeftArray[k + 1] += (mpi.rank > 0) ? -ysize * zsize * Q : (2 * mpi.size - 1) * ysize * zsize * Q;
            nodeLeftArray[k + 5] = distr.getDistribution(nodeLeft + k + 5);
            nodeLeftArray[k + 5] += (mpi.rank > 0) ? -ysize * zsize * Q : (2 * mpi.size - 1) * ysize * zsize * Q;
            nodeLeftArray[k + 7] = distr.getDistribution(nodeLeft + k + 7);
            nodeLeftArray[k + 7] += (mpi.rank > 0) ? -ysize * zsize * Q : (2 * mpi.size - 1) * ysize * zsize * Q;
            nodeRightArray[k + 2] = distr.getDistribution(nodeRight + k + 2);
            nodeRightArray[k + 2] +=
                (mpi.rank < mpi.size - 1) ? ysize * zsize * Q : (1 - 2 * mpi.size) * ysize * zsize * Q;
            nodeRightArray[k + 6] = distr.getDistribution(nodeRight + k + 6);
            nodeRightArray[k + 6] +=
                (mpi.rank < mpi.size - 1) ? ysize * zsize * Q : (1 - 2 * mpi.size) * ysize * zsize * Q;
            nodeRightArray[k + 8] = distr.getDistribution(nodeRight + k + 8);
            nodeRightArray[k + 8] +=
                (mpi.rank < mpi.size - 1) ? ysize * zsize * Q : (1 - 2 * mpi.size) * ysize * zsize * Q;
        }

    // 3. DISTRIBUITION COMMUNICATION
    lattice.communicateDistribution(distr);

    // 4. RETRIEVING DATA FROM DISTRIBUTION AFTER COMMUNICATION
    std::vector<double> distrNodeLeft(distr.getDistributionPointer(width * ysize * zsize),
                                      distr.getDistributionPointer((width + 1) * ysize * zsize));
    std::vector<double> distrNodeRight(distr.getDistributionPointer((xExt - width - 1) * ysize * zsize),
                                       distr.getDistributionPointer((xExt - width) * ysize * zsize));

    // 5. COMPARISON WITH THE EXPECTED VALUES
    EXPECT_TRUE(ArraysMatch(distrNodeLeft, nodeLeftArray));
    EXPECT_TRUE(ArraysMatch(distrNodeRight, nodeRightArray));
}

//
// D3Q15
//
TEST(ParallelX_Test, TestCommunicateDistributionD3Q15) {
    Lattice::ResetParallelTracking();

    Lattice lattice;
    using DQ = D3Q15;
    constexpr int Q = DQ::Q;
    constexpr int neighbors = Parallel_Pattern::mNumDirections * 2;
    DataOldNew<Lattice, DQ> data;
    int xExt = Lattice::LXdiv;  // length of local buffer including external halo region along X
    int nodeLeft = (width - 1) * ysize * zsize * Q;
    int nodeRight = (xExt - width) * ysize * zsize * Q;

    auto distr = data.getDistributionObject();
    ASSERT_EQ(ysize, Lattice::LYdiv);
    ASSERT_EQ(zsize, Lattice::LZdiv);
    ASSERT_EQ(distr.mv_Distribution.size(), xExt * ysize * zsize * Q);  // #elements in the distribution
    std::fill(distr.mv_Distribution.begin(), distr.mv_Distribution.end(), 0.0);

    // 1. DISTRIBUITION INITIALISATION FOR COMMUNICATION
    for (int i = 0; i < neighbors; ++i)
        for (int y = 0; y < ysize; ++y)  // #elements along Y axis
            for (int z = 0; z < zsize; ++z)
                for (int j = 0; j < Q; ++j)  // #directions in the stencil
                {
                    int k = (i == 0)
                                ? nodeLeft + (z + y * zsize) * Q + j
                                : nodeRight + (z + y * zsize) * Q + j;  // element in the distribution (with directions)
                    double value = ((neighbors * mpi.rank + i) * ysize * zsize + (z + y * zsize)) * Q +
                                   j;  // the value in the element of the distribution
                    distr.getDistribution(k) = value;
                }

    // 2. SETTING THE EXPECTED VECTOR WITH EXPECTED VALUES AFTER COMMUNICATION
    std::vector<double> nodeLeftArray(Q * ysize * zsize, 0.0);
    std::vector<double> nodeRightArray(Q * ysize * zsize, 0.0);
    for (int y = 0; y < ysize; ++y)
        for (int z = 0; z < zsize; ++z) {
            int k = (z + y * zsize) * Q;
            nodeLeftArray[k + 1] = distr.getDistribution(nodeLeft + k + 1);
            nodeLeftArray[k + 1] += (mpi.rank > 0) ? -ysize * zsize * Q : (2 * mpi.size - 1) * ysize * zsize * Q;
            nodeLeftArray[k + 7] = distr.getDistribution(nodeLeft + k + 7);
            nodeLeftArray[k + 7] += (mpi.rank > 0) ? -ysize * zsize * Q : (2 * mpi.size - 1) * ysize * zsize * Q;
            nodeLeftArray[k + 9] = distr.getDistribution(nodeLeft + k + 9);
            nodeLeftArray[k + 9] += (mpi.rank > 0) ? -ysize * zsize * Q : (2 * mpi.size - 1) * ysize * zsize * Q;
            nodeLeftArray[k + 11] = distr.getDistribution(nodeLeft + k + 11);
            nodeLeftArray[k + 11] += (mpi.rank > 0) ? -ysize * zsize * Q : (2 * mpi.size - 1) * ysize * zsize * Q;
            nodeLeftArray[k + 14] = distr.getDistribution(nodeLeft + k + 14);
            nodeLeftArray[k + 14] += (mpi.rank > 0) ? -ysize * zsize * Q : (2 * mpi.size - 1) * ysize * zsize * Q;
            nodeRightArray[k + 2] = distr.getDistribution(nodeRight + k + 2);
            nodeRightArray[k + 2] +=
                (mpi.rank < mpi.size - 1) ? ysize * zsize * Q : (1 - 2 * mpi.size) * ysize * zsize * Q;
            nodeRightArray[k + 8] = distr.getDistribution(nodeRight + k + 8);
            nodeRightArray[k + 8] +=
                (mpi.rank < mpi.size - 1) ? ysize * zsize * Q : (1 - 2 * mpi.size) * ysize * zsize * Q;
            nodeRightArray[k + 10] = distr.getDistribution(nodeRight + k + 10);
            nodeRightArray[k + 10] +=
                (mpi.rank < mpi.size - 1) ? ysize * zsize * Q : (1 - 2 * mpi.size) * ysize * zsize * Q;
            nodeRightArray[k + 12] = distr.getDistribution(nodeRight + k + 12);
            nodeRightArray[k + 12] +=
                (mpi.rank < mpi.size - 1) ? ysize * zsize * Q : (1 - 2 * mpi.size) * ysize * zsize * Q;
            nodeRightArray[k + 13] = distr.getDistribution(nodeRight + k + 13);
            nodeRightArray[k + 13] +=
                (mpi.rank < mpi.size - 1) ? ysize * zsize * Q : (1 - 2 * mpi.size) * ysize * zsize * Q;
        }

    // 3. DISTRIBUITION COMMUNICATION
    lattice.communicateDistribution(distr);

    // 4. RETRIEVING DATA FROM DISTRIBUTION AFTER COMMUNICATION
    std::vector<double> distrNodeLeft(distr.getDistributionPointer(width * ysize * zsize),
                                      distr.getDistributionPointer((width + 1) * ysize * zsize));
    std::vector<double> distrNodeRight(distr.getDistributionPointer((xExt - width - 1) * ysize * zsize),
                                       distr.getDistributionPointer((xExt - width) * ysize * zsize));

    // 5. COMPARISON WITH THE EXPECTED VALUES
    EXPECT_TRUE(ArraysMatch(distrNodeLeft, nodeLeftArray));
    EXPECT_TRUE(ArraysMatch(distrNodeRight, nodeRightArray));
}
