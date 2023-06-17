#include "test_main.hh"
#include "Collide.hh"
#include "Forcing.hh"

#include "Lattice.hh"
#include "Data.hh"
#include "Parallel.hh"
#include "Stencil.hh"

double epsilon = 1e-10;

using Lattice = LatticeProperties<Data1, NoParallel, 1, 1>;

TEST(CollideTest, computeGammaD2Q9) {
  CollisionBase<Lattice,D2Q9> collision;
  double v[2] = {1,2};
  double vv = v[0]*v[0] + v[1]*v[1];
  for (int i=0; i<9; i++) {
    int c[2] = {D2Q9::Ci_x[i], D2Q9::Ci_y[i]};
    double vc = v[0]*c[0] + v[1]*c[1];
    double gamma = D2Q9::Weights[i] * (1 + 3*vc + 4.5*vc*vc - 1.5*vv);
    EXPECT_NEAR(collision.computeGamma(v, i), gamma, epsilon);
  }
}


TEST(CollideTest, computeZerothMomentD2Q9) {
  CollisionBase<Lattice,D2Q9> collision;
  double distr[9] = {1,1,1,1,1,1,1,1,1};
  EXPECT_NEAR(collision.computeZerothMoment(distr), 9, epsilon);
}


TEST(CollideTest, computeFirstMomentD2Q9) {
  CollisionBase<Lattice,D2Q9> collision;
  double distr[9] = {1,1,0,1,0,1,0,1,0};
  EXPECT_NEAR(collision.computeFirstMoment(distr,0), 3, epsilon);
  EXPECT_NEAR(collision.computeFirstMoment(distr,1), 1, epsilon);
}


TEST(CollideTest, collideSRTD2Q9) {
  CollisionBase<Lattice,D2Q9> collision;
  double old = 2;
  double eq = 1;
  EXPECT_NEAR(collision.collide<SRT>(old, eq, 1), eq, epsilon);
  EXPECT_NEAR(collision.collide<SRT>(old, eq, 0.5), 0.5*(eq+old), epsilon);
  EXPECT_NEAR(collision.collide<SRT>(old, eq, 0.9), 0.9*eq+0.1*old, epsilon);
}

TEST(CollideTest, forceGuoSRTD2Q9) {
  Guo<Lattice,D2Q9> force;
  Velocity<Lattice> vel;  
  InverseTau<Lattice> itau;
  //CollisionBase<Lattice,D2Q9> collision;
  double tau = 0.9;
  itau.getParameter(0)=1./tau;
  double f[2] = {1, 2};
  std::copy(f, f+2, force.ma_Force);
  double v[2] = {-1, 1};
  std::copy(v, v+2, &vel.getParameter()[0]);
  double vf = v[0]*f[0] + v[1]*f[1];
  for (int i=0; i<9; i++) {
    double factor = (1 - 1/(2*tau)) * D2Q9::Weights[i];
    int c[2] = {D2Q9::Ci_x[i], D2Q9::Ci_y[i]};
    double cf = c[0]*f[0] + c[1]*f[1];
    double cv = c[0]*v[0] + c[1]*v[1];
    double fi = factor * (3*cf + 9*cv*cf - 3*vf);
    EXPECT_NEAR(force.compute(i,0), fi, epsilon);
  }
}
