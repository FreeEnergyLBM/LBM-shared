#include "test_main.hh"
#include "Lattice.hh"
// #include "Data.hh"
#include "Parallel.hh"

using Lattice = LatticeProperties<Data1, X_Parallel<1>, 6, 1>;


TEST(ParallelTest, TestCommunicate) {
  Lattice lattice;
  Density<Lattice> density;
  ASSERT_EQ(density.mv_Parameter.size(), 5);
  if (mpi.rank == 0) {
    density.mv_Parameter = {0, 1, 2, 3, 0};
  } else if (mpi.rank == 1) {
    density.mv_Parameter = {0, 4, 5, 6, 0};
  }

  lattice.communicate(density);
  std::vector<double> densityArray = (mpi.rank==0) ? std::vector<double>{6,1,2,3,4} : std::vector<double>{3,4,5,6,1};
  EXPECT_TRUE(ArraysMatch(density.mv_Parameter, densityArray));
}


TEST(ParallelTest, TestCommunicateDistribution) {
  Lattice lattice;
  Data1<Lattice,D2Q5> data;

  auto distr = data.getDistributionObject();
  ASSERT_EQ(distr.mv_Distribution.size(), 5*5);
  for (int i=0; i<2; i++) {
    for (int j=0; j<5; j++) {
      int k = (i==0) ? j : 4*5+j;
      double value = 5*(2*mpi.rank + i) + j;
      distr.mv_Distribution[k] = value;
    }
  }
  lattice.communicateDistribution(distr);

  std::vector<double> distrNode2(distr.getDistributionPointer(1), distr.getDistributionPointer(2));
  std::vector<double> distrNode4(distr.getDistributionPointer(3), distr.getDistributionPointer(4));
  std::vector<double> node2Array = (mpi.rank==0) ? std::vector<double>{0,16,0,0,0} : std::vector<double>{0,6,0,0,0};
  std::vector<double> node4Array = (mpi.rank==0) ? std::vector<double>{0,0,12,0,0} : std::vector<double>{0,0,2,0,0};
  EXPECT_TRUE(ArraysMatch(distrNode2, node2Array));
  EXPECT_TRUE(ArraysMatch(distrNode4, node4Array));
}
