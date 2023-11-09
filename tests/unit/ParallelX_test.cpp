// XPROCS 2
#include "test_main.hh"
#include "Lattice.hh"
#include "Parallel.hh"

constexpr int xsize = 8;
constexpr int ysize = 8;
constexpr int width = 2;
using Parallel_Pattern = ParallelX<width>;
using Lattice = LatticeProperties<Parallel_Pattern, xsize, ysize>;


TEST(ParallelX_Test, TestCommunicate) {
  Lattice lattice;
  std::vector<double> &density = Density<>::getInstance<Lattice>().mv_Parameter;
  int xLocalSize = Lattice::LXdiv;
  int yLocalSize = Lattice::LYdiv;
  ASSERT_EQ(density.size(), xLocalSize*yLocalSize);
  std::fill(density.begin(), density.end(), 0.0);
  for(int x=width; x<xLocalSize-width; ++x)
    for(int y=0; y<yLocalSize; ++y) {
      int i = y + x*yLocalSize;
      density[i] = (mpi.rank*(xLocalSize-2*width)-width)*ysize + i + 1;
    }
  std::vector<double> densityArray(density.begin(), density.end());
  for(int x=0; x<width; ++x)
    for(int y=0; y<yLocalSize; ++y) {
      int i = y + x*yLocalSize;
      densityArray[i] = density[i+width*yLocalSize];
      densityArray[i] += (mpi.rank>0) ? -width*yLocalSize : xsize*ysize-width*yLocalSize;
    }
  for(int x=xLocalSize-width; x<xLocalSize; ++x)
    for(int y=0; y<yLocalSize; ++y) {
      int i = y + x*yLocalSize;
      densityArray[i] = density[i-width*yLocalSize];
      densityArray[i] += (mpi.rank<mpi.size-1) ? width*yLocalSize : width*yLocalSize-xsize*ysize;
    }

  lattice.communicate(Density<>::getInstance<Lattice>());
  EXPECT_TRUE(ArraysMatch(density, densityArray));
}


TEST(ParallelX_Test, TestCommunicateDistributionD2Q5) {
  Lattice lattice;
  using DQ = D2Q5;
  constexpr int Q = DQ::Q;
  constexpr int neighbors = Parallel_Pattern::mNumDirections*2;
  DataOldNew<Lattice,DQ> data;
  int xLocalSize = Lattice::LXdiv; // length of local buffer including external halo region along X
  int yLocalSize = Lattice::LYdiv;
  int nodeLeft = (width-1)*yLocalSize*Q;
  int nodeRight = (xLocalSize-width)*yLocalSize*Q;

  auto distr = data.getDistributionObject();
  ASSERT_EQ(distr.mv_Distribution.size(), xLocalSize*yLocalSize*Q); // #elements in the distribution
  std::fill(distr.mv_Distribution.begin(), distr.mv_Distribution.end(), 0.0);
  for (int i=0; i<neighbors; i++) {
    for (int y=0; y<yLocalSize; y++) // #elements along Y axis
      for (int j=0; j<Q; j++) { // #directions in the stencil
        int k = (i==0) ? nodeLeft+y*Q+j : nodeRight+y*Q+j; // element idx in the distribution (with directions)
        double value = ((neighbors*mpi.rank+i)*yLocalSize+y)*Q+j; // the value in the element of the distribution
        distr.getDistribution(k) = value;
      }
  }
  std::vector<double> nodeLeftArray(Q*yLocalSize, 0.0);
  std::vector<double> nodeRightArray(Q*yLocalSize, 0.0);
  for(int y=0; y<yLocalSize; ++y){
    int k = y*Q;
    nodeLeftArray[k+1]   = distr.getDistribution(nodeLeft+k+1);
    nodeLeftArray[k+1]  += (mpi.rank>0)          ? -yLocalSize*Q : (4*ysize-yLocalSize)*Q;
    nodeRightArray[k+2]  = distr.getDistribution(nodeRight+k+2);
    nodeRightArray[k+2] += (mpi.rank<mpi.size-1) ?  yLocalSize*Q : (yLocalSize-4*ysize)*Q;
  }

  lattice.communicateDistribution(distr);

  std::vector<double> distrNodeLeft(distr.getDistributionPointer(width*yLocalSize), distr.getDistributionPointer((width+1)*yLocalSize));
  std::vector<double> distrNodeRight(distr.getDistributionPointer((xLocalSize-width-1)*yLocalSize), distr.getDistributionPointer((xLocalSize-width)*yLocalSize));

  EXPECT_TRUE(ArraysMatch(distrNodeLeft, nodeLeftArray));
  EXPECT_TRUE(ArraysMatch(distrNodeRight, nodeRightArray));
}


TEST(ParallelX_Test, TestCommunicateDistributionD2Q9) {
  Lattice lattice;
  using DQ = D2Q9;
  constexpr int Q = DQ::Q;
  constexpr int neighbors = Parallel_Pattern::mNumDirections*2;
  DataOldNew<Lattice,DQ> data;
  int xLocalSize = Lattice::LXdiv; // length of local buffer including external halo region along X
  int yLocalSize = Lattice::LYdiv;
  int faceSize = yLocalSize * Lattice::LZdiv;

  auto distr = data.getDistributionObject();
  ASSERT_EQ(distr.mv_Distribution.size(), xLocalSize*yLocalSize*Q); // #elements in the distribution
  std::fill(distr.mv_Distribution.begin(), distr.mv_Distribution.end(), 0.0);
  for (int i=0; i<neighbors; i++) {
    for (int y=0; y<yLocalSize; y++) // #elements along Y axis
      for (int j=0; j<Q; j++) { // #directions in the stencil
        int k = (i==0) ? ((width-1)*faceSize+y)*Q+j : ((xLocalSize-width)*faceSize+y)*Q+j; // element idx in the distribution (with directions)
        double value = ((neighbors*mpi.rank+i)*yLocalSize+y)*Q+j; // the value in the element of the distribution
        distr.mv_Distribution[k] = value;
      }
  }
  std::vector<double> nodeLeftArray(Q*yLocalSize, 0.0);
  std::vector<double> nodeRightArray(Q*yLocalSize, 0.0);
  for(int y=0; y<yLocalSize; ++y){
    nodeLeftArray[y*Q+1] = distr.mv_Distribution[((width-1)*faceSize+y)*Q+1];
    nodeLeftArray[y*Q+1] += (mpi.rank>0) ? -yLocalSize*Q : (4*ysize-yLocalSize)*Q;
    nodeLeftArray[y*Q+5] = distr.mv_Distribution[((width-1)*faceSize+y)*Q+5];
    nodeLeftArray[y*Q+5] += (mpi.rank>0) ? -yLocalSize*Q : (4*ysize-yLocalSize)*Q;
    nodeLeftArray[y*Q+7] = distr.mv_Distribution[((width-1)*faceSize+y)*Q+7];
    nodeLeftArray[y*Q+7] += (mpi.rank>0) ? -yLocalSize*Q : (4*ysize-yLocalSize)*Q;
    nodeRightArray[y*Q+2] = distr.mv_Distribution[((xLocalSize-width)*faceSize+y)*Q+2];
    nodeRightArray[y*Q+2] += (mpi.rank<mpi.size-1) ? yLocalSize*Q : (yLocalSize-4*ysize)*Q;
    nodeRightArray[y*Q+6] = distr.mv_Distribution[((xLocalSize-width)*faceSize+y)*Q+6];
    nodeRightArray[y*Q+6] += (mpi.rank<mpi.size-1) ? yLocalSize*Q : (yLocalSize-4*ysize)*Q;
    nodeRightArray[y*Q+8] = distr.mv_Distribution[((xLocalSize-width)*faceSize+y)*Q+8];
    nodeRightArray[y*Q+8] += (mpi.rank<mpi.size-1) ? yLocalSize*Q : (yLocalSize-4*ysize)*Q;
  }

  lattice.communicateDistribution(distr);

  std::vector<double> distrNodeLeft(distr.getDistributionPointer(width*yLocalSize), distr.getDistributionPointer((width+1)*yLocalSize));
  std::vector<double> distrNodeRight(distr.getDistributionPointer((xLocalSize-width-1)*yLocalSize), distr.getDistributionPointer((xLocalSize-width)*yLocalSize));
  EXPECT_TRUE(ArraysMatch(distrNodeLeft, nodeLeftArray));
  EXPECT_TRUE(ArraysMatch(distrNodeRight, nodeRightArray));
}
