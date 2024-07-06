#include "AddOns/Gradients.hh"

#include "Parameters.hh"
#include "Trait.hh"
#include "test_main.hh"

using Lattice = LatticeProperties<NoParallel, 3, 1>;
using Trait = DefaultTrait<Lattice>;

TEST(GradientsMulti, result) {
    auto orderParam1 = OrderParameter<0>::getAddress<Lattice>(0);
    auto orderParam2 = OrderParameter<1>::getAddress<Lattice>(0);
    orderParam1[0] = -1;
    orderParam1[1] = -1;
    orderParam1[2] = 1;
    orderParam2[0] = 1;
    orderParam2[1] = -1;
    orderParam2[2] = 1;

    GradientsMulti<std::tuple<OrderParameter<0>, OrderParameter<1>>, std::tuple<CentralXYZ, LaplacianCentral>>
        gradients;
    gradients.template compute<Trait>(1);

    // Check gradient
    double gradOP1 = Gradient<OrderParameter<0>>::template get<Lattice>(1);
    double gradOP2 = Gradient<OrderParameter<1>>::template get<Lattice>(1);
    EXPECT_FLOAT_EQ(gradOP1, 1);
    EXPECT_FLOAT_EQ(gradOP2, 0);

    // Check laplacian
    double lapOP1 = Laplacian<OrderParameter<0>>::template get<Lattice>(1);
    double lapOP2 = Laplacian<OrderParameter<1>>::template get<Lattice>(1);
    EXPECT_FLOAT_EQ(lapOP1, 2);
    EXPECT_FLOAT_EQ(lapOP2, 4);
}

TEST(GradientsMultiInstance, result) {
    auto orderParam1 = OrderParameter<0>::getAddress<Lattice>(0);
    auto orderParam2 = OrderParameter<1>::getAddress<Lattice>(0);
    orderParam1[0] = -1;
    orderParam1[1] = -1;
    orderParam1[2] = 1;
    orderParam2[0] = 1;
    orderParam2[1] = -1;
    orderParam2[2] = 1;

    GradientsMultiInstance<OrderParameter, 2, CentralXYZ, LaplacianCentral> gradients;
    gradients.template compute<Trait>(1);

    // Check gradient
    double gradOP1 = Gradient<OrderParameter<0>>::template get<Lattice>(1);
    double gradOP2 = Gradient<OrderParameter<1>>::template get<Lattice>(1);
    EXPECT_FLOAT_EQ(gradOP1, 1);
    EXPECT_FLOAT_EQ(gradOP2, 0);

    // Check laplacian
    double lapOP1 = Laplacian<OrderParameter<0>>::template get<Lattice>(1);
    double lapOP2 = Laplacian<OrderParameter<1>>::template get<Lattice>(1);
    EXPECT_FLOAT_EQ(lapOP1, 2);
    EXPECT_FLOAT_EQ(lapOP2, 4);
}

TEST(GradientsMultiParam, result) {
    auto orderParam1 = OrderParameter<0>::getAddress<Lattice>(0);
    auto orderParam2 = OrderParameter<1>::getAddress<Lattice>(0);
    orderParam1[0] = -1;
    orderParam1[1] = -1;
    orderParam1[2] = 1;
    orderParam2[0] = 1;
    orderParam2[1] = -1;
    orderParam2[2] = 1;

    GradientsMultiParam<LaplacianCentral, OrderParameter<0>, OrderParameter<1>> gradients;
    gradients.template compute<Trait>(1);

    // Check laplacian
    double lapOP1 = Laplacian<OrderParameter<0>>::template get<Lattice>(1);
    double lapOP2 = Laplacian<OrderParameter<1>>::template get<Lattice>(1);
    EXPECT_FLOAT_EQ(lapOP1, 2);
    EXPECT_FLOAT_EQ(lapOP2, 4);
}

TEST(GradientsMultiStencil, result) {
    auto orderParam = OrderParameter<0>::getAddress<Lattice>(0);
    orderParam[0] = -1;
    orderParam[1] = -1;
    orderParam[2] = 1;

    GradientsMultiStencil<OrderParameter<0>, CentralXYZ, LaplacianCentral> gradients;
    gradients.template compute<Trait>(1);

    // Check
    double gradOP = Gradient<OrderParameter<0>>::template get<Lattice>(1);
    double lapOP = Laplacian<OrderParameter<0>>::template get<Lattice>(1);
    EXPECT_FLOAT_EQ(gradOP, 1);
    EXPECT_FLOAT_EQ(lapOP, 2);
}
