#include "LBModels/Binary.hh"

#include "test_main.hh"

TEST(Binary, computeConcentration) {
    using Lattice = LatticeProperties<NoParallel, 3>;
    Binary<Lattice> model;

    OrderParameter<>::template initialise<Lattice>(-1, 0);
    EXPECT_EQ(model.computeConcentration(0, 0), 0);
    EXPECT_EQ(model.computeConcentration(0, 1), 1);

    OrderParameter<>::template initialise<Lattice>(0, 1);
    EXPECT_EQ(model.computeConcentration(1, 0), 0.5);
    EXPECT_EQ(model.computeConcentration(1, 1), 0.5);

    OrderParameter<>::template initialise<Lattice>(1, 2);
    EXPECT_EQ(model.computeConcentration(2, 0), 1);
    EXPECT_EQ(model.computeConcentration(2, 1), 0);
}

TEST(FlowFieldBinary, computeConcentration) {
    using Lattice = LatticeProperties<NoParallel, 3>;
    FlowFieldBinary<Lattice> model;

    OrderParameter<>::template initialise<Lattice>(-1, 0);
    EXPECT_EQ(model.computeConcentration(0, 0), 0);
    EXPECT_EQ(model.computeConcentration(0, 1), 1);

    OrderParameter<>::template initialise<Lattice>(0, 1);
    EXPECT_EQ(model.computeConcentration(1, 0), 0.5);
    EXPECT_EQ(model.computeConcentration(1, 1), 0.5);

    OrderParameter<>::template initialise<Lattice>(1, 2);
    EXPECT_EQ(model.computeConcentration(2, 0), 1);
    EXPECT_EQ(model.computeConcentration(2, 1), 0);
}

