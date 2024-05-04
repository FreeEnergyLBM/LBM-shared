#include "BoundaryModels/BounceBack.hh"

#include "Data.hh"
#include "Global.hh"
#include "LBModels/ModelBase.hh"
#include "Lattice.hh"
#include "Parallel.hh"
#include "test_main.hh"

using Lattice = LatticeProperties<NoParallel, 2, 1>;

TEST(BounceBackTest, TestNodePair) {
    using Trait = DefaultTrait<Lattice>::SetStencil<D2Q9>;

    BoundaryLabels<Lattice::NDIM>::get<Lattice>(0).Id = 1;
    BounceBack bb;
    bb.setNodeID(1);

    DataOldNew<Lattice, D2Q9> data;
    auto& distr = data.getDistributionObject();
    distr.mv_Distribution = {0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 0, 0, 0, 0};

    bb.compute<Trait>(distr, 0);
    std::vector<double> distrNode1(distr.getDistributionPointer(1), distr.getDistributionPointer(2));
    EXPECT_TRUE(ArraysMatch(distrNode1, {0, 2, 1, 0, 0, 6, 5, 8, 7}));
}
