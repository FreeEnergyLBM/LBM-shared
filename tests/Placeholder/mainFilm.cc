#include <math.h>
#include <stdlib.h>

#include <lbm.hh>

// This script simulates a four component layered poiseuille flow setup.
// You can modify the densities and relaxation times via the 'setDensities' and 'setTaus' functions.
// You can modify the body force magnitude in the setMagnitudeX function

const int lx = 2;    // Size of domain in x direction
const int ly = 100;  // Size of domain in y direction

const int timesteps = 300000;    // Number of iterations to perform
const int saveInterval = 25000;  // Interval to save global data

// Parameters to control the surface tension and width of the diffuse interface
// Use these if you want the surface tensions to all be the same
// double BETA=0.001;
// double GAMMA=-BETA*6.25;

double RADIUS = 25.0;

double A = 0.0025;
double kappa = 0.005;
const double Hsat = 0.5;

using Lattice = LatticeProperties<NoParallel, lx, ly>;
double offsetx = 0.0;
double offsety = 0.0;
// Function used to define the solid geometry
// Here we set a solid at the top and bottom, in the conditions that return 1;
int initBoundary(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1) / 2. - offsetx) * (xx - (lx - 1) / 2. - offsetx) +
                 (yy - (ly - 1) / 2. - offsety) * (yy - (ly - 1) / 2. - offsety);

    // if (yy > 1 && (yy >= ly - 1 || xx == 1 || xx == lx - 2)) return 4;
    if (yy >= ly - 2) return 4;
    if (yy <= 2) return 1;
    // if (yy == 1 || yy == ly - 2 || xx == 1 || xx == lx - 2) return 4;
    // if (xx <= 1) return 1;
    // else if (xx == lx - 2) return 4;
    // else if (xx >= lx - 1) return 1;
    // if(sqrt(rr2)<RADIUS) return 5;
    // if(sqrt(rr2)<RADIUS) return 5;
    if (yy < ly / 2. + offsety) return 5;
    return 0;
}

double initFluid(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1) / 2. - offsetx) * (xx - (lx - 1) / 2. - offsetx) +
                 (yy - (ly - 1) / 2. - offsety) * (yy - (ly - 1) / 2. - offsety);
    // return 0.25*((double)rand()/(double)RAND_MAX);

    if (yy <= 1) return 0;

    // else return 0.5-0.5*tanh(2*(sqrt(rr2)-RADIUS)/(sqrt(8*kappa/A)));
    // return 0.5-0.5*tanh(2*((xx - lx/2.-offsety))/(sqrt(8*kappa/A)));
    return 0.5 - 0.5 * tanh(2 * ((yy - ly / 2. - offsety)) / (sqrt(8 * kappa / A)));
    // if(sqrt(rr2)<RADIUS) return 1;
    // else return 0;
}

double initHumidity(int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1) / 2. - offsetx) * (xx - (lx - 1) / 2. - offsetx) +
                 (yy - (ly - 1) / 2. - offsety) * (yy - (ly - 1) / 2. - offsety);
    // return 0.25*((double)rand()/(double)RAND_MAX);
    // if(sqrt(rr2)<RADIUS&&yy > 1) return Hsat;
    // else return 0;
    if (yy <= ly / 2. + offsety) return Hsat;
    // if (xx <= 0 || xx >= lx - 1) return 0;
    // return Hsat - Hsat*(xx - (lx/2.+offsety))/(lx-lx/2.-offsety-1.5);
    return Hsat - Hsat * (yy - (ly / 2. + offsety)) / (ly - ly / 2. - offsety - 1.5);

    // else return 0;
    // int yy = computeY(ly, 1, k);
    // return 0.5*tanh(2*(yy-3*ly/4)/(D))-0*0.5*tanh(2*(yy-ly)/(D))+0.5-0*0.5*tanh(2*(yy)/(D));
}
/*
bool interfaceCondition(const double& val, int k){
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    //return 0.25*((double)rand()/(double)RAND_MAX);
    if(sqrt(rr2)<RADIUS) return false;
    else return true;
    //return val<0.5;
}
*/

bool interfaceCondition(const double& val, int k) { return val < 0.5; }

bool interfaceConditionK(int k) { return OrderParameter<>::get<Lattice>(k) < 0.5; }

// Function used to define the fluid
// Here we set a tanh transition in the y direction. This should match the equilibrium profile.
/*
double initFluid1(int k) {
    //int xx = computeXGlobal<Lattice>(k);
    //double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    //return 0.25*((double)rand()/(double)RAND_MAX);
    //return 0.5+0.5*tanh((sqrt(rr2)-RADIUS)/(sqrt(2*kappa/A)));
    int yy = computeY(ly, 1, k);
    return 0.5*tanh(2*(yy-3*ly/4)/(D))-0*0.5*tanh(2*(yy-ly)/(D))+0.5-0*0.5*tanh(2*(yy)/(D));
}

double initFluid2(int k) {
    int yy = computeY(ly, 1, k);
    return 0.5*tanh(2*(yy-ly/2)/(D))-0.5*tanh(2*(yy-3*ly/4)/(D));
}

double initFluid3(int k) {
    int yy = computeY(ly, 1, k);
    return 0.5*tanh(2*(yy-ly/4)/(D))-0.5*tanh(2*(yy-ly/2)/(D));
}
*/

using traithumid = DefaultTraitHumidity<Lattice>::SetStencil<D2Q5>;
using traitpressure =
    typename DefaultTraitPressureLeeHumidity<Lattice>::template SetBoundary<PressureOutflow, BounceBack>;
// using traitpressure = typename DefaultTraitPressureLee<Lattice> :: AddForce<BodyForce<>>;

double distancefunc(int k, int idx) {
    using Stencil = typename traithumid::Stencil;

    double normaldist = -0.5 * sqrt(8 * kappa / A) * atanh((OrderParameter<>::get<Lattice>(k) - 0.5) / 0.5);
    double* gradorderparam = GradientOrderParameter<>::getAddress<Lattice, Lattice::NDIM>(k, 0);
    // float_or_double magnitudegradient = sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1],
    // 2)+pow(gradorderparam[2], 2));
    std::vector<double> normal;
    if constexpr (Lattice::NDIM == 1)
        normal = {-gradorderparam[0] / sqrt(pow(gradorderparam[0], 2))};
    else if constexpr (Lattice::NDIM == 2)
        normal = {-gradorderparam[0] / sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2)),
                  -gradorderparam[1] / sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2))};
    else
        normal = {-gradorderparam[0] /
                      sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2) + pow(gradorderparam[2], 2)),
                  -gradorderparam[1] /
                      sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2) + pow(gradorderparam[2], 2)),
                  -gradorderparam[2] /
                      sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2) + pow(gradorderparam[2], 2))};
    double normdotci = 0;
    double magci = 0;
    for (int xyz = 0; xyz < Lattice::NDIM; xyz++) {
        normdotci += normal[xyz] * Stencil::Ci_xyz(xyz)[idx];
        magci += Stencil::Ci_xyz(xyz)[idx] * Stencil::Ci_xyz(xyz)[idx];
    }
    // return 0.5; //!!!!!!!!!!!!!!!!
    double dist = fabs(normaldist / normdotci / sqrt(magci));
    // std::cout<<dist<<" "<<normaldist<<" "<<cos(normdotci/sqrt(magci))<<" "<<normdotci<<" "<<magci<<std::endl;
    // if(dist == 0.6) std::cout<<k<<std::endl;
    // std::cout<<dist<<std::endl;
    if (idx == 0 || std::isnan(dist) || std::isinf(dist))
        return 0.5;
    else if (dist >= 1)
        return 1;
    else if (dist <= 0)
        return 0;
    return dist;  // normdotci/magci/magci;
}

int main(int argc, char** argv) {
    // mpi.init();

    // Set up the lattice, including the resolution and data/parallelisation method

    // We need to modify the traits of the navier stokes model to include a bodyforce and change the collision model to
    // MRT, so that high viscosity contrasts produce accurate results. We need to add a chemical potential calculator to
    // one of the models too, given to the first NComponent model.

    EvaporationHumidity<Lattice, traithumid> humidity;
    PressureLeeHumidity<Lattice, traitpressure> pressure;
    BinaryLeeHumidity<Lattice> binary;

    // PressureLee<Lattice,traitpressure> pressure;
    // BinaryLee<Lattice> binary;

    // BinaryLee<Lattice> binary;

    binary.setDensity1(1);
    binary.setDensity2(1);

    binary.getForce<MuSourceLocal>().setBeta(A);
    binary.getForce<MuSourceNonLocal>().setBeta(A);

    binary.getPostProcessor<ChemicalPotentialCalculatorBinaryLee>().setA(A);
    binary.getPostProcessor<ChemicalPotentialCalculatorBinaryLee>().setKappa(kappa);

    binary.getPreProcessor<MassLossCalculatorInterpolated>().setInterfaceHumidity(Hsat);
    binary.getPreProcessor<MassLossCalculatorInterpolated>().setDiffusivity(0.008);
    binary.getPreProcessor<MassLossCalculatorInterpolated>().setInterfaceWidth(sqrt(8 * kappa / A));
    binary.getPreProcessor<MassLossCalculatorInterpolated>().setPhiGasLiquid(0, 1);

    double theta = 2.0 * M_PI / 4.0;
    double wettingprefactor = -cos(theta) * sqrt(2 * A / kappa);

    // binary.getPreProcessor<GradientsWettingMultiStencil<OrderParameter<>, CentralXYZWetting, CentralQWetting,
    // MixedXYZWetting, MixedQWetting, LaplacianCentralWetting>>().setPrefactor(wettingprefactor);

    binary
        .getPreProcessor<GradientsWettingMultiStencil<OrderParameter<>, CentralXYZBounceBack, CentralQBounceBack,
                                                      MixedXYZBounceBack, MixedQBounceBack, LaplacianCentralWetting>>()
        .setPrefactor(wettingprefactor);

    // binary.getPreProcessor<MassLossCalculator>().setInterfaceHumidity(0.5);
    // binary.getPreProcessor<MassLossCalculator>().setDiffusivity(0.2);

    using dbtype = InterpolatedDirichlet;

    // humidity.getBoundary<InterpolatedDirichlet,0>().setInterfaceDistanceFunction(distancefunc);
    // humidity.getBoundary<dbtype,0>().setNodeID(5);
    // humidity.getBoundary<dbtype,0>().setInterfaceVal(0.8);

    humidity.getBoundary<InterpolatedDirichlet>().setInterfaceDistanceFunction(distancefunc);
    humidity.getBoundary<dbtype>().setNodeID(5);
    humidity.getBoundary<dbtype>().setInterfaceVal(Hsat);

    humidity.getBoundary<Dirichlet>().setNodeID(4);
    humidity.getBoundary<Dirichlet>().setInterfaceVal(0.0);

    humidity.setDiffusivity(0.008);

    humidity.getPreProcessor<HumidityBoundaryLabels>().setInterfaceCondition(interfaceCondition);
    humidity.getPreProcessor<SetHumidityLiquid>().setInterfaceVal(Hsat);

    humidity.getPostProcessor<GradientsInterface<Humidity<>, CentralXYZInterfaceBounceBack>>().setInterfaceDistance(
        distancefunc);
    // humidity.getPostProcessor<GradientsInterface<Humidity<>,
    // CentralXYZInterfaceBounceBack>>().setInterfaceCondition(interfaceConditionK);
    humidity.getPostProcessor<GradientsInterface<Humidity<>, CentralXYZInterfaceBounceBack>>().setInterfaceVal(Hsat);

    humidity.getForce<EvaporationHumiditySource<EvaporationSourceMethod>>().setInterfaceHumidity(Hsat);

    pressure.getForce<EvaporationPressureSource<EvaporationSourceMethod>>().setInterfaceHumidity(Hsat);
    // pressure.getForce<BodyForce<>>().setMagnitudeX(0.0000001);
    SaveHandler<Lattice> saver("data/");
    saver.saveHeader(timesteps, saveInterval);  // Create a header with lattice information (lx, ly, lz, NDIM (2D or
                                                // 3D), timesteps, saveInterval)

    // Define the solid and fluid using the functions above
    Geometry<Lattice>::initialiseBoundaries(initBoundary);
    OrderParameter<>::set<Lattice>(initFluid);
    OrderParameterOld<>::set<Lattice>(initFluid);
    Humidity<>::set<Lattice>(initHumidity);

    // Algorithm creates an object that can run our chosen LBM model
    // Algorithm lbm(PressureNavierStokes,NCompAllenCahn1,NCompAllenCahn2,NCompAllenCahn3);

    // Set up the handler object for saving data

    Algorithm lbm(binary, humidity);
    // Algorithm lbm(binary);
    // Algorithm lbm(pressure,binary,humidity);
    // Algorithm lbm(humidity);
    // Algorithm lbm(pressure,binary);

    // Perform the main LBM loop
    saver.saveBoundaries(0);
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        // Save the desired parameters, producing a binary file for each.
        if (timestep % saveInterval == 0) {
            if (mpi.rank == 0) std::cout << "Saving at timestep " << timestep << "." << std::endl;
            // saver.saveParameter<BoundaryLabels<>>(timestep);

            // saver.saveBoundaries(timestep);
            saver.saveParameter<Humidity<>>(timestep);
            saver.saveParameter<ChemicalPotential<>>(timestep);
            saver.saveParameter<Density<>>(timestep);
            saver.saveParameter<Pressure<>>(timestep);
            saver.saveParameter<OrderParameter<>>(timestep);
            saver.saveParameter<MassSink<>>(timestep);
            saver.saveParameter<Velocity<>, Lattice::NDIM>(timestep);
            saver.saveParameter<VelocityOld<>, Lattice::NDIM>(timestep);
            saver.saveParameter<GradientHumidity<>, Lattice::NDIM>(timestep);
        }

        // Evolve by one timestep
        lbm.evolve();
    }
    return 0;
}
