#include <math.h>
#include <stdlib.h>

#include <lbm.hh>

// This script simulates a four component layered poiseuille flow setup.
// You can modify the densities and relaxation times via the 'setDensities' and 'setTaus' functions.
// You can modify the body force magnitude in the setMagnitudeX function

int lx = 400;  // Size of domain in x direction
int ly = 100;  // Size of domain in y direction
int lz = 1;
int timesteps = 5000;     // Number of iterations to perform
int saveInterval = 1000;  // Interval to save global data
double dens1 = 1;
double dens2 = 1;
double dens3 = 1;
double radius = 30;
std::string datadir = "data/";
double s01 = 0.001;
double s02 = 0.001;
double s12 = 0.001;
double gamma0 = s01 + s02 - s12;
double gamma1 = s01 + s12 - s02;
double gamma2 = s02 + s12 - s01;
double interfacewidth = 5;
double tau1 = 1.0;
double tau2 = 1.0;
double tau3 = 1.0;
std::vector<double> gammas = {};
InputParameters params;

#ifdef MPIPARALLEL
using Lattice = LatticePropertiesRuntime<ParallelX<2>, 2>;
#else
using Lattice = LatticePropertiesRuntime<NoParallel, 2>;
#endif

double initFluid1(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, lz, k);
    double rr2 = (xx - (lx - 1) / 2.) * (xx - (lx - 1) / 2.) + (yy - (ly - 1) / 2.) * (yy - (ly - 1) / 2.);
    if (rr2 < radius * radius) {
        return 1;
    } else
        return 0;
}

double initFluid2(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1) / 2.) * (xx - (lx - 1) / 2.) + (yy - (ly - 1) / 2.) * (yy - (ly - 1) / 2.);
    if (yy >= (ly - 1) / 2 && rr2 >= radius * radius)
        return 1;
    else
        return 0;
}

void initParams(std::string inputfile) {
    params.addParameter<int>(lx, "lx");
    params.addParameter<int>(ly, "ly");
    params.addParameter<int>(lz, "lz");
    params.addParameter<int>(timesteps, "timesteps");
    params.addParameter<int>(saveInterval, "saveInterval");
    params.addParameter<double>(radius, "radius");
    params.addParameter<double>(s01, "s01");
    params.addParameter<double>(s02, "s02");
    params.addParameter<double>(s12, "s12");
    params.addParameter<double>(dens1, "dens1");
    params.addParameter<double>(dens2, "dens2");
    params.addParameter<double>(dens3, "dens3");
    params.addParameter<double>(radius, "radius");
    params.addParameter<std::string>(datadir, "datadir");
    params.addParameter<double>(tau1, "tau1");
    params.addParameter<double>(tau2, "tau2");
    params.addParameter<double>(tau3, "tau3");
    params.addParameter<double>(interfacewidth, "interfacewidth");

    params.readInput(inputfile);

    gamma0 = s01 + s02 - s12;
    gamma1 = s01 + s12 - s02;
    gamma2 = s02 + s12 - s01;
    gammas = {gamma0, gamma1, gamma2};
    Lattice::init(lx, ly, lz);
}

template <typename TTrait = DefaultTraitPressureTernaryLee<Lattice>>
auto initPressure() {
    PressureTernaryLee<Lattice> pressure;

    pressure.setCollideID({0});

    pressure.setTau1(tau1);
    pressure.setTau2(tau2);
    pressure.setTau3(tau3);

    pressure.setDensity1(dens1);
    pressure.setDensity2(dens2);
    pressure.setDensity3(dens3);

    return pressure;
}

template <int id = 0, class TTrait = DefaultTraitTernaryLee<Lattice, id>>
auto initTernary() {
    using Trait = typename std::conditional<
        id == 0, typename TTrait::template AddProcessorIdx<0, ChemicalPotentialCalculatorTernaryLee>, TTrait>::type;

    TernaryLee<Lattice, id, Trait> ternary;
    ternary.setCollideID({0});

    ternary.template getForce<MuTernarySourceLocal<id>>().setBetaInterfaceWidth(gamma0, interfacewidth);
    ternary.template getForce<MuTernarySourceNonLocal<id>>().setBetaInterfaceWidth(gamma0, interfacewidth);

    if constexpr (id == 0) {
        ternary.template getProcessor<ChemicalPotentialCalculatorTernaryLee>().setSurfaceTension(s01, s02, s12);
        ternary.template getProcessor<ChemicalPotentialCalculatorTernaryLee>().setInterfaceWidth(interfacewidth);
    }

    return ternary;
}