#include "LBModels/TernaryLee.hh"
#include "test_main.hh"

TEST(TernaryLee, computeConcentration) {
    using Lattice = LatticeProperties<NoParallel, 2>;
    TernaryLee<Lattice, 0> model;

    OrderParameter<0>::template initialise<Lattice>(0, 0);
    OrderParameter<1>::template initialise<Lattice>(1, 0);
    EXPECT_EQ(model.computeConcentration(0, 0), 0);
    EXPECT_EQ(model.computeConcentration(0, 1), 1);
    EXPECT_EQ(model.computeConcentration(0, 2), 0);

    OrderParameter<0>::template initialise<Lattice>(0, 1);
    OrderParameter<1>::template initialise<Lattice>(0, 1);
    EXPECT_EQ(model.computeConcentration(1, 0), 0);
    EXPECT_EQ(model.computeConcentration(1, 1), 0);
    EXPECT_EQ(model.computeConcentration(1, 2), 1);
}

TEST(PressureTernaryLee, computeConcentration) {
    using Lattice = LatticeProperties<NoParallel, 3>;
    PressureTernaryLee<Lattice> model;

    OrderParameter<0>::template initialise<Lattice>(0, 0);
    OrderParameter<1>::template initialise<Lattice>(1, 0);
    EXPECT_EQ(model.computeConcentration(0, 0), 0);
    EXPECT_EQ(model.computeConcentration(0, 1), 1);
    EXPECT_EQ(model.computeConcentration(0, 2), 0);

    OrderParameter<0>::template initialise<Lattice>(0, 1);
    OrderParameter<1>::template initialise<Lattice>(0, 1);
    EXPECT_EQ(model.computeConcentration(1, 0), 0);
    EXPECT_EQ(model.computeConcentration(1, 1), 0);
    EXPECT_EQ(model.computeConcentration(1, 2), 1);
}
