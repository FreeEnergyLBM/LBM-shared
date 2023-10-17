#include "test_main.hh"
#include "Global.hh"
#include "BoundaryModels/BounceBack.hh"
#include "Lattice.hh"
#include "Data.hh"
#include "Parallel.hh"
#include "LBModels/ModelBase.hh"

using Lattice = LatticeProperties<DataOldNew, NoParallel, 2, 1>;

TEST(BounceBackTest, TestNodePair) {
  using Trait = DefaultTrait<Lattice> ::SetStencil<D2Q9>;

  BoundaryLabels<Lattice::NDIM>::get<Lattice>(0).Id = 1;

  DataOldNew<Lattice,D2Q9> data;
  auto& distr = data.getDistributionObject();
  distr.mv_Distribution = {0,1,2,3,4,5,6,7,8, 0,0,0,0,0,0,0,0};

  BounceBack bb;
  bb.compute<Trait>(distr, 0);
  std::vector<double> distrNode1(distr.getDistributionPointer(1), distr.getDistributionPointer(2));
  EXPECT_TRUE(ArraysMatch(distrNode1, {0,2,1,0,0,6,5,8,7}));
}
