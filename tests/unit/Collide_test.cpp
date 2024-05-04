#include "Collide.hh"

#include "Data.hh"
#include "Forcing.hh"
#include "LBModels/ModelBase.hh"
#include "Lattice.hh"
#include "Parallel.hh"
#include "Stencil.hh"
#include "test_main.hh"

double epsilon = 1e-10;

using Lattice = LatticeProperties<NoParallel, 1, 1>;

TEST(CollideTest, computeGammaD2Q9) {
    CollisionBase<Lattice, D2Q9> collision;
    double v[2] = {1, 2};
    double vv = v[0] * v[0] + v[1] * v[1];
    for (int i = 0; i < 9; i++) {
        int c[2] = {D2Q9::Ci_x[i], D2Q9::Ci_y[i]};
        double vc = v[0] * c[0] + v[1] * c[1];
        double gamma = D2Q9::Weights[i] * (1 + 3 * vc + 4.5 * vc * vc - 1.5 * vv);
        EXPECT_NEAR(collision.computeGamma(v, i), gamma, epsilon);
    }
}

TEST(CollideTest, computeZerothMomentD2Q9) {
    CollisionBase<Lattice, D2Q9> collision;
    double distr[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    EXPECT_NEAR(collision.computeZerothMoment(distr), 9, epsilon);
}

TEST(CollideTest, computeFirstMomentD2Q9) {
    CollisionBase<Lattice, D2Q9> collision;
    double distr[9] = {1, 1, 0, 1, 0, 1, 0, 1, 0};
    EXPECT_NEAR(collision.computeFirstMoment(distr, 0), 3, epsilon);
    EXPECT_NEAR(collision.computeFirstMoment(distr, 1), 1, epsilon);
}

TEST(CollideTest, collideSRTD2Q9) {
    using Collision = SRT<D2Q9>;
    double old[1] = {2};
    double eq[1] = {1};
    EXPECT_NEAR(Collision::collide<Lattice>(old, eq, 1, 0), eq[0], epsilon);
    EXPECT_NEAR(Collision::collide<Lattice>(old, eq, 0.5, 0), 0.5 * (eq[0] + old[0]), epsilon);
    EXPECT_NEAR(Collision::collide<Lattice>(old, eq, 0.9, 0), 0.9 * eq[0] + 0.1 * old[0], epsilon);
}

// TEST(CollideTest, forceGuoSRTD2Q9) {
//   using Trait = DefaultTrait<Lattice>;
//   Guo force;
//   auto vel = Velocity<>::get<Lattice>(0);
//   auto itau = InverseTau<>::get<Lattice>(0);
//   double tau = 0.9;
//   itau = 1./tau;
//   double f[2] = {1, 2};
//   std::copy(f, f+2, force.ma_Force);
//   double v[2] = {-1, 1};
//   std::copy(v, v+2, vel);
//   double vf = v[0]*f[0] + v[1]*f[1];
//   for (int i=0; i<9; i++) {
//     double factor = (1 - 1/(2*tau)) * D2Q9::Weights[i];
//     int c[2] = {D2Q9::Ci_x[i], D2Q9::Ci_y[i]};
//     double cf = c[0]*f[0] + c[1]*f[1];
//     double cv = c[0]*v[0] + c[1]*v[1];
//     double fi = factor * (3*cf + 9*cv*cf - 3*vf);
//     EXPECT_NEAR(force.compute<Trait>(i,0), fi, epsilon);
//   }
// }
