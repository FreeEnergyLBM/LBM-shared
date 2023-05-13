#include "test_main.hh"
#include "LBModels/FlowField.hh"

#include "Lattice.hh"
#include "Data.hh"
#include "Parallel.hh"
#include "Stencil.hh"

TEST(FlowFieldTest, initialiseD2Q9) {
  FlowField<LatticeProperties<Data1,X_Parallel,1,1>> ff;
  ff.initialise();

  // Add two neighbours due to the parallel implementation
  std::vector<double> equilibrium(3*9);
  std::copy(D2Q9::Weights, D2Q9::Weights+9, equilibrium.begin()+9);

  EXPECT_TRUE(ArraysMatch(ff.getDistribution(), equilibrium));
}
