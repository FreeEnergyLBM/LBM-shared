// YPROCS 4
#include "Lattice.hh"
#include "Parallel.hh"
#include "test_main.hh"

constexpr int xsize = 128;  // global==local since no parallelisation along X
constexpr int YSIZE = 128;  // global size along Y, parallelised
constexpr int zsize = 128;  // global==local since no parallelisation along Z
constexpr int width = 2;
using Parallel_Pattern = ParallelY<width>;
using Lattice = LatticeProperties<Parallel_Pattern, xsize, YSIZE, zsize>;

// In the tests below:
// - I introduce the yExt variable ('y extended')
//   since parallelisation is done along Y and the local size along Y includes 'width' from both sides.
// - 'xsize' and 'zsize' remain the same full sizes since parallelisation isn't done along X and Z.

// In the 'TestCommunicateParameter' TEST:
// - I construct the 2 vectors:
//   (1) for communication ('density') and
//   (2) for testing against (1) ('densityArray' containing the expected values).
// - The communication vectors are constructed such that its elements have consecutive numbers beginning from 1
//   (ordered firstly along Z, secondly along Y, thirdly along X)

TEST(ParallelY_Test, TestCommunicateParameter) {
    Lattice::ResetParallelTracking();

    Lattice lattice;
    std::vector<double> &density = Density<>::getInstance<Lattice>().mv_Parameter;
    int yExt = Lattice::LYdiv;
    ASSERT_EQ(xsize, Lattice::LXdiv);
    ASSERT_EQ(zsize, Lattice::LZdiv);
    ASSERT_EQ(density.size(), xsize * yExt * zsize);
    std::fill(density.begin(), density.end(), 0.0);
    for (int x = 0; x < xsize; ++x)
        for (int y = width; y < yExt - width; ++y)
            for (int z = 0; z < zsize; ++z) {
                int i = z + y * zsize + x * yExt * zsize;
                // rank=0:  z + (y-width)*zsize + x*YSIZE*zsize + 1
                // rank>0:  z + (y-width)*zsize + x*YSIZE*zsize + mpi.rank*(yExt-2*width)*zsize + 1
                //     or:  z + ((y-width) + mpi.rank*(yExt-2*width) + x*YSIZE)*zsize + 1
                density[i] = i + x * (YSIZE - yExt) * zsize + (mpi.rank * (yExt - 2 * width) - width) * zsize + 1;
                //      z + y*zsize + x*YSIZE*zsize + (mpi.rank*(yExt-2*width)-width)*zsize + 1
                // or:  z + y*zsize + x*yExt*zsize + x*(YSIZE-yExt)*zsize + (mpi.rank*(yExt-2*width)-width)*zsize + 1
            }
    std::vector<double> densityArray(density.begin(), density.end());
    for (int x = 0; x < xsize; ++x)
        for (int y = 0; y < width; ++y)
            for (int z = 0; z < zsize; ++z) {
                int i = z + y * zsize + x * yExt * zsize;
                densityArray[i] = density[i + width * zsize];
                densityArray[i] += (mpi.rank > 0) ? -width * zsize : (YSIZE - width) * zsize;
            }
    for (int x = 0; x < xsize; ++x)
        for (int y = yExt - width; y < yExt; ++y)
            for (int z = 0; z < zsize; ++z) {
                int i = z + y * zsize + x * yExt * zsize;
                densityArray[i] = density[i - width * zsize];
                densityArray[i] += (mpi.rank < mpi.size - 1) ? width * zsize : (width - YSIZE) * zsize;
            }

    lattice.communicate(Density<>::getInstance<Lattice>());

    EXPECT_TRUE(ArraysMatch(density, densityArray));
}

// In the 'TestCommunicateDistribution*' TESTs:
// - I construct the 2 sets of vectors:
//   (1) for communication ('distrNodeBottom' and 'distrNodeTop') and
//   (2) for testing against (1) ('nodeBottomArray' and 'nodeTopArray' containing the expected values).
// - The communication vectors are constructed such that its elements have consecutive numbers beginning from 0
//   (ordered firstly along velocity directions, and then along Z, along Y, and along X)

//
// D2Q5
//
TEST(ParallelY_Test, TestCommunicateDistributionD2Q5) {
    Lattice::ResetParallelTracking();

    Lattice lattice;
    using DQ = D2Q5;
    constexpr int Q = DQ::Q;
    constexpr int neighbors = Parallel_Pattern::mNumDirections * 2;
    DataOldNew<Lattice, DQ> data;
    int yExt = Lattice::LYdiv;
    int nodeBottom = (width - 1) * zsize * Q;  // + (z+x*yExt*zsize)
    int nodeTop = (yExt - width) * zsize * Q;

    auto distr = data.getDistributionObject();
    ASSERT_EQ(xsize, Lattice::LXdiv);
    ASSERT_EQ(zsize, Lattice::LZdiv);
    ASSERT_EQ(distr.mv_Distribution.size(), xsize * yExt * zsize * Q);  // #elements in the distribution
    std::fill(distr.mv_Distribution.begin(), distr.mv_Distribution.end(), 0.0);

    // 1. DISTRIBUITION INITIALISATION FOR COMMUNICATION
    for (int i = 0; i < neighbors; ++i) {
        for (int x = 0; x < xsize; ++x)
            for (int z = 0; z < zsize; ++z) {
                int l = (z + x * yExt * zsize) * Q;
                for (int j = 0; j < Q; ++j)  // #directions in the stencil
                {
                    int k = (i == 0) ? nodeBottom + l + j
                                     : nodeTop + l + j;  // element k in the distribution (with directions)
                    double value = ((neighbors * mpi.rank + i) * xsize * zsize + (z + x * zsize)) * Q +
                                   j;  // the value in the element of the distribution
                    distr.getDistribution(k) = value;
                }
            }
    }

    // 2. SETTING THE EXPECTED VECTOR WITH EXPECTED VALUES AFTER COMMUNICATION
    std::vector<double> nodeBottomArray(Q * xsize * zsize, 0.0);
    std::vector<double> nodeTopArray(Q * xsize * zsize, 0.0);
    for (int x = 0; x < xsize; ++x)
        for (int z = 0; z < zsize; ++z) {
            int k = (z + x * zsize) * Q;
            int l = (z + x * yExt * zsize) * Q;
            nodeBottomArray[k + 3] = distr.getDistribution(nodeBottom + l + 3);
            nodeBottomArray[k + 3] += (mpi.rank > 0) ? -xsize * zsize * Q : (2 * mpi.size - 1) * xsize * zsize * Q;
            nodeTopArray[k + 4] = distr.getDistribution(nodeTop + l + 4);
            nodeTopArray[k + 4] +=
                (mpi.rank < mpi.size - 1) ? xsize * zsize * Q : (1 - 2 * mpi.size) * xsize * zsize * Q;
        }

    // 3. DISTRIBUITION COMMUNICATION
    lattice.communicateDistribution(distr);

    // 4. RETRIEVING DATA FROM DISTRIBUTION AFTER COMMUNICATION
    std::vector<double> distrNodeBottom(xsize * zsize * Q, 0.0);
    std::vector<double> distrNodeTop(xsize * zsize * Q, 0.0);
    for (int x = 0; x < xsize; ++x)
        for (int z = 0; z < zsize; ++z) {
            int k = (z + x * zsize) * Q;  //(x*yExt + width)*zsize*Q
            int l = (z + x * yExt * zsize) * Q;
            for (int j = 0; j < Q; ++j) {  // #directions in the stencil
                distrNodeBottom[k + j] = distr.getDistribution(nodeBottom + l + Q * zsize + j);
                distrNodeTop[k + j] = distr.getDistribution(nodeTop + l - Q * zsize + j);
            }
        }

    // 5. COMPARISON WITH THE EXPECTED VALUES
    EXPECT_TRUE(ArraysMatch(distrNodeBottom, nodeBottomArray));
    EXPECT_TRUE(ArraysMatch(distrNodeTop, nodeTopArray));
}

//
// D2Q9
//
TEST(ParallelY_Test, TestCommunicateDistributionD2Q9) {
    Lattice::ResetParallelTracking();

    Lattice lattice;
    using DQ = D2Q9;
    constexpr int Q = DQ::Q;
    constexpr int neighbors = Parallel_Pattern::mNumDirections * 2;
    DataOldNew<Lattice, DQ> data;
    int yExt = Lattice::LYdiv;
    int nodeBottom = (width - 1) * zsize * Q;  // + (z+x*yExt*zsize)
    int nodeTop = (yExt - width) * zsize * Q;

    auto distr = data.getDistributionObject();
    ASSERT_EQ(xsize, Lattice::LXdiv);
    ASSERT_EQ(zsize, Lattice::LZdiv);
    ASSERT_EQ(distr.mv_Distribution.size(), xsize * yExt * zsize * Q);  // #elements in the distribution
    std::fill(distr.mv_Distribution.begin(), distr.mv_Distribution.end(), 0.0);

    // 1. DISTRIBUITION INITIALISATION FOR COMMUNICATION
    for (int i = 0; i < neighbors; ++i)
        for (int x = 0; x < xsize; ++x)  // #elements along X axis
            for (int z = 0; z < zsize; ++z) {
                int l = (z + x * yExt * zsize) * Q;
                for (int j = 0; j < Q; ++j)  // #directions in the stencil
                {
                    int k = (i == 0) ? nodeBottom + l + j
                                     : nodeTop + l + j;  // element k in the distribution (with directions)
                    double value = ((neighbors * mpi.rank + i) * xsize * zsize + (z + x * zsize)) * Q +
                                   j;  // the value in the element of the distribution
                    distr.getDistribution(k) = value;
                }
            }

    // 2. SETTING THE EXPECTED VECTOR WITH EXPECTED VALUES AFTER COMMUNICATION
    std::vector<double> nodeBottomArray(Q * xsize * zsize, 0.0);
    std::vector<double> nodeTopArray(Q * xsize * zsize, 0.0);
    for (int x = 0; x < xsize; ++x)
        for (int z = 0; z < zsize; ++z) {
            int k = (z + x * zsize) * Q;
            int l = (z + x * yExt * zsize) * Q;
            nodeBottomArray[k + 3] = distr.getDistribution(nodeBottom + l + 3);
            nodeBottomArray[k + 3] += (mpi.rank > 0) ? -xsize * zsize * Q : (2 * mpi.size - 1) * xsize * zsize * Q;
            nodeBottomArray[k + 5] = distr.getDistribution(nodeBottom + l + 5);
            nodeBottomArray[k + 5] += (mpi.rank > 0) ? -xsize * zsize * Q : (2 * mpi.size - 1) * xsize * zsize * Q;
            nodeBottomArray[k + 8] = distr.getDistribution(nodeBottom + l + 8);
            nodeBottomArray[k + 8] += (mpi.rank > 0) ? -xsize * zsize * Q : (2 * mpi.size - 1) * xsize * zsize * Q;
            nodeTopArray[k + 4] = distr.getDistribution(nodeTop + l + 4);
            nodeTopArray[k + 4] +=
                (mpi.rank < mpi.size - 1) ? xsize * zsize * Q : (1 - 2 * mpi.size) * xsize * zsize * Q;
            nodeTopArray[k + 6] = distr.getDistribution(nodeTop + l + 6);
            nodeTopArray[k + 6] +=
                (mpi.rank < mpi.size - 1) ? xsize * zsize * Q : (1 - 2 * mpi.size) * xsize * zsize * Q;
            nodeTopArray[k + 7] = distr.getDistribution(nodeTop + l + 7);
            nodeTopArray[k + 7] +=
                (mpi.rank < mpi.size - 1) ? xsize * zsize * Q : (1 - 2 * mpi.size) * xsize * zsize * Q;
        }

    // 3. DISTRIBUITION COMMUNICATION
    lattice.communicateDistribution(distr);

    // 4. RETRIEVING DATA FROM DISTRIBUTION AFTER COMMUNICATION
    std::vector<double> distrNodeBottom(xsize * zsize * Q, 0.0);
    std::vector<double> distrNodeTop(xsize * zsize * Q, 0.0);
    for (int x = 0; x < xsize; ++x)
        for (int z = 0; z < zsize; ++z) {
            int k = (z + x * zsize) * Q;  //(x*yExt + width)*zsize*Q
            int l = (z + x * yExt * zsize) * Q;
            for (int j = 0; j < Q; ++j)  // #directions in the stencil
            {
                distrNodeBottom[k + j] = distr.getDistribution(nodeBottom + l + Q * zsize + j);
                distrNodeTop[k + j] = distr.getDistribution(nodeTop + l - Q * zsize + j);
            }
        }

    // 5. COMPARISON WITH THE EXPECTED VALUES
    EXPECT_TRUE(ArraysMatch(distrNodeBottom, nodeBottomArray));
    EXPECT_TRUE(ArraysMatch(distrNodeTop, nodeTopArray));
}

//
// D3Q15
//
TEST(ParallelY_Test, TestCommunicateDistributionD3Q15) {
    Lattice::ResetParallelTracking();

    Lattice lattice;
    using DQ = D3Q15;
    constexpr int Q = DQ::Q;
    constexpr int neighbors = Parallel_Pattern::mNumDirections * 2;
    DataOldNew<Lattice, DQ> data;
    int yExt = Lattice::LYdiv;
    int nodeBottom = (width - 1) * zsize * Q;  // + (z+x*yExt*zsize)
    int nodeTop = (yExt - width) * zsize * Q;

    auto distr = data.getDistributionObject();
    ASSERT_EQ(xsize, Lattice::LXdiv);
    ASSERT_EQ(zsize, Lattice::LZdiv);
    ASSERT_EQ(distr.mv_Distribution.size(), xsize * yExt * zsize * Q);  // #elements in the distribution
    std::fill(distr.mv_Distribution.begin(), distr.mv_Distribution.end(), 0.0);

    // 1. DISTRIBUITION INITIALISATION FOR COMMUNICATION
    for (int i = 0; i < neighbors; ++i)
        for (int x = 0; x < xsize; ++x)  // #elements along X axis
            for (int z = 0; z < zsize; ++z) {
                int l = (z + x * yExt * zsize) * Q;
                for (int j = 0; j < Q; ++j)  // #directions in the stencil
                {
                    int k = (i == 0) ? nodeBottom + l + j
                                     : nodeTop + l + j;  // element k in the distribution (with directions)
                    double value = ((neighbors * mpi.rank + i) * xsize * zsize + (z + x * zsize)) * Q +
                                   j;  // the value in the element of the distribution
                    distr.getDistribution(k) = value;
                }
            }

    // 2. SETTING THE EXPECTED VECTOR WITH EXPECTED VALUES AFTER COMMUNICATION
    std::vector<double> nodeBottomArray(Q * xsize * zsize, 0.0);
    std::vector<double> nodeTopArray(Q * xsize * zsize, 0.0);
    for (int x = 0; x < xsize; ++x)
        for (int z = 0; z < zsize; ++z) {
            int k = (z + x * zsize) * Q;
            int l = (z + x * yExt * zsize) * Q;
            nodeBottomArray[k + 3] = distr.getDistribution(nodeBottom + l + 3);
            nodeBottomArray[k + 3] += (mpi.rank > 0) ? -xsize * zsize * Q : (2 * mpi.size - 1) * xsize * zsize * Q;
            nodeBottomArray[k + 7] = distr.getDistribution(nodeBottom + l + 7);
            nodeBottomArray[k + 7] += (mpi.rank > 0) ? -xsize * zsize * Q : (2 * mpi.size - 1) * xsize * zsize * Q;
            nodeBottomArray[k + 9] = distr.getDistribution(nodeBottom + l + 9);
            nodeBottomArray[k + 9] += (mpi.rank > 0) ? -xsize * zsize * Q : (2 * mpi.size - 1) * xsize * zsize * Q;
            nodeBottomArray[k + 12] = distr.getDistribution(nodeBottom + l + 12);
            nodeBottomArray[k + 12] += (mpi.rank > 0) ? -xsize * zsize * Q : (2 * mpi.size - 1) * xsize * zsize * Q;
            nodeBottomArray[k + 13] = distr.getDistribution(nodeBottom + l + 13);
            nodeBottomArray[k + 13] += (mpi.rank > 0) ? -xsize * zsize * Q : (2 * mpi.size - 1) * xsize * zsize * Q;
            nodeTopArray[k + 4] = distr.getDistribution(nodeTop + l + 4);
            nodeTopArray[k + 4] +=
                (mpi.rank < mpi.size - 1) ? xsize * zsize * Q : (1 - 2 * mpi.size) * xsize * zsize * Q;
            nodeTopArray[k + 8] = distr.getDistribution(nodeTop + l + 8);
            nodeTopArray[k + 8] +=
                (mpi.rank < mpi.size - 1) ? xsize * zsize * Q : (1 - 2 * mpi.size) * xsize * zsize * Q;
            nodeTopArray[k + 10] = distr.getDistribution(nodeTop + l + 10);
            nodeTopArray[k + 10] +=
                (mpi.rank < mpi.size - 1) ? xsize * zsize * Q : (1 - 2 * mpi.size) * xsize * zsize * Q;
            nodeTopArray[k + 11] = distr.getDistribution(nodeTop + l + 11);
            nodeTopArray[k + 11] +=
                (mpi.rank < mpi.size - 1) ? xsize * zsize * Q : (1 - 2 * mpi.size) * xsize * zsize * Q;
            nodeTopArray[k + 14] = distr.getDistribution(nodeTop + l + 14);
            nodeTopArray[k + 14] +=
                (mpi.rank < mpi.size - 1) ? xsize * zsize * Q : (1 - 2 * mpi.size) * xsize * zsize * Q;
        }

    // 3. DISTRIBUITION COMMUNICATION
    lattice.communicateDistribution(distr);

    // 4. RETRIEVING DATA FROM DISTRIBUTION AFTER COMMUNICATION
    std::vector<double> distrNodeBottom(xsize * zsize * Q, 0.0);
    std::vector<double> distrNodeTop(xsize * zsize * Q, 0.0);
    for (int x = 0; x < xsize; ++x)
        for (int z = 0; z < zsize; ++z) {
            int k = (z + x * zsize) * Q;  //(x*yExt + width)*zsize*Q
            int l = (z + x * yExt * zsize) * Q;
            for (int j = 0; j < Q; ++j)  // #directions in the stencil
            {
                distrNodeBottom[k + j] = distr.getDistribution(nodeBottom + l + Q * zsize + j);
                distrNodeTop[k + j] = distr.getDistribution(nodeTop + l - Q * zsize + j);
            }
        }

    // 5. COMPARISON WITH THE EXPECTED VALUES
    EXPECT_TRUE(ArraysMatch(distrNodeBottom, nodeBottomArray));
    EXPECT_TRUE(ArraysMatch(distrNodeTop, nodeTopArray));
}
