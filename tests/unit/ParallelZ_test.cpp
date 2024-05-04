// ZPROCS 4
#include "Lattice.hh"
#include "Parallel.hh"
#include "test_main.hh"

constexpr int xsize = 128;  // global==local since no parallelisation along X
constexpr int ysize = 128;  // global==local since no parallelisation along Y
constexpr int ZSIZE = 128;  // global size along Z, parallelised
constexpr int width = 2;
using Parallel_Pattern = ParallelZ<width>;
using Lattice = LatticeProperties<Parallel_Pattern, xsize, ysize, ZSIZE>;

// In the tests below:
// - I introduce the zExt variable ('z extended')
//   since parallelisation is done along Z and the local size along Z includes 'width' from both sides.
// - 'xsize' and 'ysize' remain the same full sizes since parallelisation isn't done along X and Y.

// In the 'TestCommunicateParameter' TEST:
// - I construct the 2 vectors:
//   (1) for communication ('density') and
//   (2) for testing against (1) ('densityArray' containing the expected values).
// - The communication vectors are constructed such that its elements have consecutive numbers beginning from 1
//   (ordered firstly along Z, secondly along Y, thirdly along X)

TEST(ParallelZ_Test, TestCommunicateParameter) {
    Lattice lattice;
    std::vector<double> &density = Density<>::getInstance<Lattice>().mv_Parameter;
    int zExt = Lattice::LZdiv;
    ASSERT_EQ(xsize, Lattice::LXdiv);
    ASSERT_EQ(ysize, Lattice::LYdiv);
    ASSERT_EQ(density.size(), xsize * ysize * zExt);
    std::fill(density.begin(), density.end(), 0.0);
    for (int x = 0; x < xsize; ++x)
        for (int y = 0; y < ysize; ++y)
            for (int z = width; z < zExt - width; ++z) {
                int i = z + y * zExt + x * ysize * zExt;
                // rank=0:  (z-width) + y*ZSIZE + x*ysize*ZSIZE + 1
                // rank>0:  (z-width) + y*ZSIZE + x*ysize*ZSIZE + mpi.rank*(zExt-2*width) + 1
                //     or:  (z-width) + (y+x*ysize)*ZSIZE + mpi.rank*(zExt-2*width) + 1
                density[i] =
                    i + y * (ZSIZE - zExt) + x * ysize * (ZSIZE - zExt) + (mpi.rank * (zExt - 2 * width) - width) + 1;
                //      z + y*ZSIZE + x*ysize*ZSIZE + (mpi.rank*(zExt-2*width)-width) + 1
                // or:  z + y*zExt + x*ysize*zExt + y*(ZSIZE-zExt) + x*ysize*(ZSIZE-zExt) +
                // (mpi.rank*(zExt-2*width)-width) + 1
            }
    std::vector<double> densityArray(density.begin(), density.end());
    for (int x = 0; x < xsize; ++x)
        for (int y = 0; y < ysize; ++y)
            for (int z = 0; z < width; ++z) {
                int i = z + y * zExt + x * ysize * zExt;
                densityArray[i] = density[i + width];
                densityArray[i] += (mpi.rank > 0) ? -width : ZSIZE - width;
            }
    for (int x = 0; x < xsize; ++x)
        for (int y = 0; y < ysize; ++y)
            for (int z = zExt - width; z < zExt; ++z) {
                int i = z + y * zExt + x * ysize * zExt;
                densityArray[i] = density[i - width];
                densityArray[i] += (mpi.rank < mpi.size - 1) ? width : width - ZSIZE;
            }

    lattice.communicate(Density<>::getInstance<Lattice>());

    EXPECT_TRUE(ArraysMatch(density, densityArray));
}

// In the 'TestCommunicateDistribution*' TESTs:
// - I construct the 2 sets of vectors:
//   (1) for communication ('distrNodeFront' and 'distrNodeBack') and
//   (2) for testing against (1) ('nodeFrontArray' and 'nodeBackArray' containing the expected values).
// - The communication vectors are constructed such that its elements have consecutive numbers beginning from 0
//   (ordered firstly along velocity directions, and then along Z, along Y, and along X)

// Skipping tests for stencils {D2Q5, D2Q9} as there's nothing to communicate there
// since Z directions are missing in these stencils.
// That's why I've added tests for new stencil {D3Q15} which weren't tested before yet.

//
// D3Q15
//
TEST(ParallelZ_Test, TestCommunicateDistributionD3Q15) {
    Lattice::ResetParallelTracking();

    Lattice lattice;
    using DQ = D3Q15;
    constexpr int Q = DQ::Q;
    constexpr int neighbors = Parallel_Pattern::mNumDirections * 2;
    DataOldNew<Lattice, DQ> data;
    int zExt = Lattice::LZdiv;
    int nodeFront = (width - 1) * Q;
    int nodeBack = (zExt - width) * Q;

    auto distr = data.getDistributionObject();
    ASSERT_EQ(xsize, Lattice::LXdiv);
    ASSERT_EQ(ysize, Lattice::LYdiv);
    ASSERT_EQ(distr.mv_Distribution.size(), xsize * ysize * zExt * Q);
    std::fill(distr.mv_Distribution.begin(), distr.mv_Distribution.end(), 0.0);

    // 1. DISTRIBUITION INITIALISATION FOR COMMUNICATION
    for (int i = 0; i < neighbors; ++i)
        for (int x = 0; x < xsize; ++x)
            for (int y = 0; y < ysize; ++y) {
                int l = (y + x * ysize) * zExt * Q;
                for (int j = 0; j < Q; ++j) {
                    int k = (i == 0) ? nodeFront + l + j
                                     : nodeBack + l + j;  // element k in the distribution (with directions)
                    double value = ((neighbors * mpi.rank + i) * xsize * ysize + (y + x * ysize)) * Q +
                                   j;  // the value in the element of the distribution
                    distr.getDistribution(k) = value;
                }
            }

    // 2. SETTING THE EXPECTED VECTOR WITH EXPECTED VALUES AFTER COMMUNICATION
    std::vector<double> nodeFrontArray(Q * xsize * ysize, 0.0);
    std::vector<double> nodeBackArray(Q * xsize * ysize, 0.0);
    for (int x = 0; x < xsize; ++x)
        for (int y = 0; y < ysize; ++y) {
            int k = (y + x * ysize) * Q;
            int l = (y + x * ysize) * zExt * Q;
            nodeFrontArray[k + 5] = distr.getDistribution(nodeFront + l + 5);
            nodeFrontArray[k + 5] += (mpi.rank > 0) ? -xsize * ysize * Q : (2 * mpi.size - 1) * xsize * ysize * Q;
            nodeFrontArray[k + 7] = distr.getDistribution(nodeFront + l + 7);
            nodeFrontArray[k + 7] += (mpi.rank > 0) ? -xsize * ysize * Q : (2 * mpi.size - 1) * xsize * ysize * Q;
            nodeFrontArray[k + 10] = distr.getDistribution(nodeFront + l + 10);
            nodeFrontArray[k + 10] += (mpi.rank > 0) ? -xsize * ysize * Q : (2 * mpi.size - 1) * xsize * ysize * Q;
            nodeFrontArray[k + 11] = distr.getDistribution(nodeFront + l + 11);
            nodeFrontArray[k + 11] += (mpi.rank > 0) ? -xsize * ysize * Q : (2 * mpi.size - 1) * xsize * ysize * Q;
            nodeFrontArray[k + 13] = distr.getDistribution(nodeFront + l + 13);
            nodeFrontArray[k + 13] += (mpi.rank > 0) ? -xsize * ysize * Q : (2 * mpi.size - 1) * xsize * ysize * Q;
            nodeBackArray[k + 6] = distr.getDistribution(nodeBack + l + 6);
            nodeBackArray[k + 6] +=
                (mpi.rank < mpi.size - 1) ? xsize * ysize * Q : (1 - 2 * mpi.size) * xsize * ysize * Q;
            nodeBackArray[k + 8] = distr.getDistribution(nodeBack + l + 8);
            nodeBackArray[k + 8] +=
                (mpi.rank < mpi.size - 1) ? xsize * ysize * Q : (1 - 2 * mpi.size) * xsize * ysize * Q;
            nodeBackArray[k + 9] = distr.getDistribution(nodeBack + l + 9);
            nodeBackArray[k + 9] +=
                (mpi.rank < mpi.size - 1) ? xsize * ysize * Q : (1 - 2 * mpi.size) * xsize * ysize * Q;
            nodeBackArray[k + 12] = distr.getDistribution(nodeBack + l + 12);
            nodeBackArray[k + 12] +=
                (mpi.rank < mpi.size - 1) ? xsize * ysize * Q : (1 - 2 * mpi.size) * xsize * ysize * Q;
            nodeBackArray[k + 14] = distr.getDistribution(nodeBack + l + 14);
            nodeBackArray[k + 14] +=
                (mpi.rank < mpi.size - 1) ? xsize * ysize * Q : (1 - 2 * mpi.size) * xsize * ysize * Q;
        }

    // 3. DISTRIBUITION COMMUNICATION
    lattice.communicateDistribution(distr);

    // 4. RETRIEVING DATA FROM DISTRIBUTION AFTER COMMUNICATION
    std::vector<double> distrNodeFront(xsize * ysize * Q, 0.0);
    std::vector<double> distrNodeBack(xsize * ysize * Q, 0.0);
    for (int x = 0; x < xsize; ++x)
        for (int y = 0; y < ysize; ++y) {
            int k = (y + x * ysize) * Q;
            int l = (y + x * ysize) * zExt * Q;
            for (int j = 0; j < Q; ++j)  // #directions in the stencil
            {
                distrNodeFront[k + j] = distr.getDistribution(nodeFront + l + Q + j);
                distrNodeBack[k + j] = distr.getDistribution(nodeBack + l - Q + j);
            }
        }

    // 5. COMPARISON WITH THE EXPECTED VALUES
    EXPECT_TRUE(ArraysMatch(distrNodeFront, nodeFrontArray));
    EXPECT_TRUE(ArraysMatch(distrNodeBack, nodeBackArray));
}
