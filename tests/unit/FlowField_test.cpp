#include "test_main.hh"
#include "LBModels/FlowField.hh"

#include "Lattice.hh"
#include "Data.hh"
#include "Parallel.hh"
#include "Stencil.hh"

TEST(FlowFieldTest, initialiseD2Q9) {
  FlowField<LatticeProperties<NoParallel,2,2>> ff;
  ff.initialise();

  std::vector<double> equilibrium(D2Q9::Weights, D2Q9::Weights+9);
  std::vector<double> distribution(&ff.getDistribution()[0], &ff.getDistribution()[0]+9);
  EXPECT_TRUE(ArraysMatch(distribution, equilibrium));
}
