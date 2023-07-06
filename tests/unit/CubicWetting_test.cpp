#include "test_main.hh"
#include "AddOns/CubicWetting.hh"

#include "LBModels/ModelBase.hh"

using Lattice = LatticeProperties<Data1, NoParallel, 2, 1>;
using Trait = DefaultTrait<Lattice>;

TEST(CubicWetting, orderParameter) {
  SolidLabels<>::get<Lattice>(0) = true;
  auto orderParam = OrderParameter<>::getAddress<Lattice>(0);
  orderParam[1] = 1;

  // Check default (no wetting)
  CubicWetting wetting;
  wetting.compute<Trait>(0);
  wetting.compute<Trait>(1);
  EXPECT_FLOAT_EQ(orderParam[0], 1);
  EXPECT_FLOAT_EQ(orderParam[1], 1);

  // Set a contact angle (60 degrees)
  wetting.setTheta(M_PI/3);
  wetting.compute<Trait>(0);
  EXPECT_FLOAT_EQ(orderParam[0], 1);
  orderParam[1] = 0;
  wetting.compute<Trait>(0);
  EXPECT_FLOAT_EQ(orderParam[0], 1.0/(2*sqrt(2.0)));

  // Set a contact angle in degrees
  wetting.setThetaDegrees(120);
  wetting.compute<Trait>(0);
  EXPECT_FLOAT_EQ(orderParam[0], -1.0/(2*sqrt(2.0)));
}
