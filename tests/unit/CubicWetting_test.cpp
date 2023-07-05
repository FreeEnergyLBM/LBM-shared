#include "test_main.hh"
#include "AddOns/CubicWetting.hh"

using Lattice = LatticeProperties<Data1, NoParallel, 2, 1>;

TEST(CubicWetting, orderParameter) {
  SolidLabels<Lattice> solid;
  solid.getParameter(0) = true;
  OrderParameter<Lattice> orderParam;
  orderParam.getParameter(1) = 1;

  // Check default (no wetting)
  CubicWetting<Lattice,D2Q9> wetting;
  wetting.postprocess(0);
  wetting.postprocess(1);
  EXPECT_FLOAT_EQ(orderParam.getParameter(0), 1);
  EXPECT_FLOAT_EQ(orderParam.getParameter(1), 1);

  // Set a contact angle (60 degrees)
  wetting.setTheta(M_PI/3);
  wetting.postprocess(0);
  EXPECT_FLOAT_EQ(orderParam.getParameter(0), 1);
  orderParam.getParameter(1) = 0;
  wetting.postprocess(0);
  EXPECT_FLOAT_EQ(orderParam.getParameter(0), 1.0/(2*sqrt(2.0)));

  // Set a contact angle in degrees
  wetting.setThetaDegrees(120);
  wetting.postprocess(0);
  EXPECT_FLOAT_EQ(orderParam.getParameter(0), -1.0/(2*sqrt(2.0)));
}
