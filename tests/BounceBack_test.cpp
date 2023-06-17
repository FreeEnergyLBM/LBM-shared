#include "test_main.hh"
#include "Global.hh"
#include "BoundaryModels/BounceBack.hh"
#include "Lattice.hh"
#include "Data.hh"
#include "Parallel.hh"

using Lattice = LatticeProperties<Data1, NoParallel, 2, 1>;

TEST(BounceBackTest, TestNodePair) {
  Data1<Lattice,D2Q9> data;
  auto distr = data.getDistributionObject();
  distr.mv_Distribution = {0,1,2,3,4,5,6,7,8, 0,0,0,0,0,0,0,0};

  BounceBack<Lattice> bb;
  for (int iQ=0; iQ<9; iQ++) bb.compute(distr, 0, iQ);
  std::vector<double> distrNode1(distr.getDistributionPointer(1), distr.getDistributionPointer(2));
  EXPECT_TRUE(ArraysMatch(distrNode1, {0,2,1,0,0,6,5,8,7}));
}
