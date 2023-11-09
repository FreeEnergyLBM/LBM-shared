// YPROCS 2
#include "test_main.hh"
#include "Lattice.hh"
#include "Parallel.hh"

constexpr int xsize = 8;
constexpr int ysize = 8;
constexpr int width = 2;
using Parallel_Pattern = ParallelY<width>;
using Lattice = LatticeProperties<Parallel_Pattern, xsize, ysize>;


TEST(ParallelY_Test, TestCommunicate) {
  Lattice lattice;
  std::vector<double> &density = Density<>::getInstance<Lattice>().mv_Parameter;
  int xLocalSize = Lattice::LXdiv;
  int yLocalSize = Lattice::LYdiv;
  ASSERT_EQ(density.size(), xLocalSize*yLocalSize);
  std::fill(density.begin(), density.end(), 0.0);
  for(int x=0; x<xLocalSize; ++x)
    for(int y=width; y<yLocalSize-width; ++y) {
      int i = y + x*yLocalSize;
      density[i] = x*ysize + mpi.rank*(yLocalSize-2*width) + (y-width) + 1;
    }
  std::vector<double> densityArray(density.begin(), density.end());
  for(int x=0; x<xLocalSize; ++x)
    for(int y=0; y<width; ++y) {
      int i = y + x*yLocalSize;
      densityArray[i] = density[i+width];
      densityArray[i] += (mpi.rank>0) ? -width : ysize-width;
    }
  for(int x=0; x<xLocalSize; ++x)
    for(int y=yLocalSize-width; y<yLocalSize; ++y) {
      int i = y + x*yLocalSize;
      densityArray[i] = density[i-width];
      densityArray[i] += (mpi.rank<mpi.size-1) ? width : width-ysize;
    }

  lattice.communicate(Density<>::getInstance<Lattice>());

  EXPECT_TRUE(ArraysMatch(density, densityArray));
}


TEST(ParallelY_Test, TestCommunicateDistributionD2Q5) {
  Lattice lattice;
  using DQ = D2Q5;
  constexpr int Q = DQ::Q;
  constexpr int neighbors = Parallel_Pattern::mNumDirections*2;
  DataOldNew<Lattice,DQ> data;
  int xLocalSize = Lattice::LXdiv; // length of local buffer including external halo region along X
  int yLocalSize = Lattice::LYdiv;
  int nodeBottom = (width-1)*Q;
  int nodeTop = (yLocalSize-width)*Q;

  auto distr = data.getDistributionObject();
  ASSERT_EQ(distr.mv_Distribution.size(), yLocalSize*xLocalSize*Q); // #elements in the distribution
  for (int i=0; i<neighbors; i++) {
    for (int x=0; x<xLocalSize; x++) // #elements along X axis
      for (int j=0; j<Q; j++) { // #directions in the stencil
        int k = (i==0) ? nodeBottom+x*yLocalSize*Q+j : nodeTop+x*yLocalSize*Q+j; // element k in the distribution (with directions)
        double value = ((neighbors*mpi.rank+i)*xLocalSize+x)*Q+j; // the value in the element of the distribution
        distr.getDistribution(k) = value;
      }
  }
  std::vector<double> nodeBottomArray(Q*xLocalSize, 0.0);
  std::vector<double> nodeTopArray(Q*xLocalSize, 0.0);
  for(int x=0; x<xLocalSize; ++x){
    int k = x*Q;
    nodeBottomArray[k+3] = distr.getDistribution(k*yLocalSize+nodeBottom+3);
    nodeBottomArray[k+3] += (mpi.rank>0) ? -xLocalSize*Q : (4*xsize-xLocalSize)*Q;
    nodeTopArray[k+4] = distr.getDistribution(k*yLocalSize+nodeTop+4);
    nodeTopArray[k+4] += (mpi.rank<mpi.size-1) ? xLocalSize*Q : (xLocalSize-4*xsize)*Q;
  }

  lattice.communicateDistribution(distr);

  std::vector<double> distrNodeBottom(xLocalSize*Q, 0.0);
  std::vector<double> distrNodeTop(xLocalSize*Q, 0.0);
  for(int x=0; x<xLocalSize; ++x) {
    int k = x*Q;//(x*yLocalSize + width)*Q;
    for (int j=0; j<Q; j++) { // #directions in the stencil
      distrNodeBottom[k+j] = distr.getDistribution(k*yLocalSize+nodeBottom+Q+j);
      distrNodeTop[k+j] = distr.getDistribution(k*yLocalSize+nodeTop-Q+j);
    }
  }
  EXPECT_TRUE(ArraysMatch(distrNodeBottom, nodeBottomArray));
  EXPECT_TRUE(ArraysMatch(distrNodeTop, nodeTopArray));
}


TEST(ParallelY_Test, TestCommunicateDistributionD2Q9) {
  Lattice lattice;
  using DQ = D2Q9;
  constexpr int Q = DQ::Q;
  constexpr int neighbors = Parallel_Pattern::mNumDirections*2;
  DataOldNew<Lattice,DQ> data;
  int xLocalSize = Lattice::LXdiv; // length of local buffer including external halo region along X
  int yLocalSize = Lattice::LYdiv;
  int nodeBottom = (width-1)*Q;
  int nodeTop = (yLocalSize-width)*Q;

  auto distr = data.getDistributionObject();
  ASSERT_EQ(distr.mv_Distribution.size(), yLocalSize*xLocalSize*Q); // #elements in the distribution
  for (int i=0; i<neighbors; i++) {
    for (int x=0; x<xLocalSize; x++) // #elements along X axis
      for (int j=0; j<Q; j++) { // #directions in the stencil
        int k = (i==0) ? nodeBottom+x*yLocalSize*Q+j : nodeTop+x*yLocalSize*Q+j; // element k in the distribution (with directions)
        double value = ((neighbors*mpi.rank+i)*xLocalSize+x)*Q+j; // the value in the element of the distribution
        distr.getDistribution(k) = value;
      }
  }
  std::vector<double> nodeBottomArray(Q*xLocalSize, 0.0);
  std::vector<double> nodeTopArray(Q*xLocalSize, 0.0);
  for(int x=0; x<xLocalSize; ++x){
    int k = x*Q;
    nodeBottomArray[k+3] = distr.getDistribution(k*yLocalSize+nodeBottom+3);
    nodeBottomArray[k+3] += (mpi.rank>0) ? -xLocalSize*Q : (4*xsize-xLocalSize)*Q;
    nodeBottomArray[k+5] = distr.getDistribution(k*yLocalSize+nodeBottom+5);
    nodeBottomArray[k+5] += (mpi.rank>0) ? -xLocalSize*Q : (4*xsize-xLocalSize)*Q;
    nodeBottomArray[k+8] = distr.getDistribution(k*yLocalSize+nodeBottom+8);
    nodeBottomArray[k+8] += (mpi.rank>0) ? -xLocalSize*Q : (4*xsize-xLocalSize)*Q;
    nodeTopArray[k+4] = distr.getDistribution(k*yLocalSize+nodeTop+4);
    nodeTopArray[k+4] += (mpi.rank<mpi.size-1) ? xLocalSize*Q : (xLocalSize-4*xsize)*Q;
    nodeTopArray[k+6] = distr.getDistribution(k*yLocalSize+nodeTop+6);
    nodeTopArray[k+6] += (mpi.rank<mpi.size-1) ? xLocalSize*Q : (xLocalSize-4*xsize)*Q;
    nodeTopArray[k+7] = distr.getDistribution(k*yLocalSize+nodeTop+7);
    nodeTopArray[k+7] += (mpi.rank<mpi.size-1) ? xLocalSize*Q : (xLocalSize-4*xsize)*Q;
  }

  lattice.communicateDistribution(distr);

  std::vector<double> distrNodeBottom(xLocalSize*Q, 0.0);
  std::vector<double> distrNodeTop(xLocalSize*Q, 0.0);
  for(int x=0; x<xLocalSize; ++x) {
    int k = x*Q;//(x*yLocalSize + width)*Q;
    for (int j=0; j<Q; j++) { // #directions in the stencil
      distrNodeBottom[k+j] = distr.getDistribution(k*yLocalSize+nodeBottom+Q+j);
      distrNodeTop[k+j] = distr.getDistribution(k*yLocalSize+nodeTop-Q+j);
    }
  }
  EXPECT_TRUE(ArraysMatch(distrNodeBottom, nodeBottomArray));
  EXPECT_TRUE(ArraysMatch(distrNodeTop, nodeTopArray));
}
