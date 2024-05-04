#include <math.h>
#include <stdlib.h>

#include <lbm.hh>

// This script simulates a four component layered poiseuille flow setup.
// You can modify the densities and relaxation times via the 'setDensities' and 'setTaus' functions.
// You can modify the body force magnitude in the setMagnitudeX function

const int lx = 5;    // Size of domain in x direction
const int ly = 100;  // Size of domain in y direction

const int timesteps = 200000;    // Number of iterations to perform
const int saveInterval = 10000;  // Interval to save global data

// Parameters to control the surface tension and width of the diffuse interface
// Use these if you want the surface tensions to all be the same
// double BETA=0.001;
// double GAMMA=-BETA*6.25;

double RADIUS = 20.0;

double A = 0.0025;
double kappa = A * 3.125;  // 0.005;
const double Hsat = 0.2;

using Lattice = LatticeProperties<NoParallel, lx, ly>;
double offsetx = 0.0;
double offsety = -ly / 2.;
// Function used to define the solid geometry
// Here we set a solid at the top and bottom, in the conditions that return 1;
int initBoundary(const int k) {
    // int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    // double rr2 = (xx - (lx-1)/2. - offsetx) * (xx - (lx-1)/2. - offsetx) + (yy - (ly-1)/2. - offsety) * (yy -
    // (ly-1)/2. - offsety);

    if (yy <= 0 || yy >= ly - 1) return 1;
    return 0;
}

double initFluid(const int k) {
    // int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    // double rr2 = (xx - (lx-1)/2. - offsetx) * (xx - (lx-1)/2. - offsetx) + (yy - (ly-1)/2. - offsety) * (yy -
    // (ly-1)/2. - offsety);

    // if (yy <= 0 || yy >= ly - 1) return 0;
    return 0.5 - 0.5 * tanh(2 * ((yy - ly / 2.)) / (sqrt(8 * kappa / A)));
}

using traitpressure = typename DefaultTraitPressureLee<Lattice>::AddForce<BodyForce<>>::SetCollisionOperator<MRT>;

int main(int argc, char **argv) {
    // mpi.init();

    // Set up the lattice, including the resolution and data/parallelisation method

    // We need to modify the traits of the navier stokes model to include a bodyforce and change the collision model to
    // MRT, so that high viscosity contrasts produce accurate results. We need to add a chemical potential calculator to
    // one of the models too, given to the first NComponent model.

    PressureLee<Lattice, traitpressure> pressure;
    BinaryLee<Lattice> binary;

    binary.setDensity1(1);
    binary.setDensity2(1);

    binary.getForce<MuSourceLocal>().setBeta(A);
    binary.getForce<MuSourceNonLocal>().setBeta(A);

    binary.getPostProcessor<ChemicalPotentialCalculatorBinaryLee>().setA(A);
    binary.getPostProcessor<ChemicalPotentialCalculatorBinaryLee>().setKappa(kappa);

    double theta = M_PI / 4.0;
    double wettingprefactor = -cos(theta) * sqrt(2 * A / kappa);

    binary
        .getPostProcessor<GradientsWettingMultiStencil<OrderParameter<>, CentralXYZWetting, CentralQWetting,
                                                       MixedXYZWetting, MixedQWetting, LaplacianCentralWetting>>()
        .setPrefactor(wettingprefactor);

    pressure.getForce<BodyForce<>>().setMagnitudeX(0.000000);  // 1);
    SaveHandler<Lattice> saver("data/");
    saver.saveHeader(timesteps, saveInterval);  // Create a header with lattice information (lx, ly, lz, NDIM (2D or
                                                // 3D), timesteps, saveInterval)

    // Define the solid and fluid using the functions above
    Geometry<Lattice>::initialiseBoundaries(initBoundary);
    OrderParameter<>::set<Lattice>(initFluid);

    // Algorithm creates an object that can run our chosen LBM model

    // Set up the handler object for saving data

    Algorithm lbm(binary, pressure);

    // Perform the main LBM loop
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        // Save the desired parameters, producing a binary file for each.
        if (timestep % saveInterval == 0) {
            if (mpi.rank == 0) std::cout << "Saving at timestep " << timestep << "." << std::endl;
            // saver.saveParameter<BoundaryLabels<>>(timestep);
            saver.saveBoundaries(timestep);
            saver.saveParameter<ChemicalPotential<>>(timestep);
            saver.saveParameter<Density<>>(timestep);
            saver.saveParameter<Pressure<>>(timestep);
            saver.saveParameter<OrderParameter<>>(timestep);
            saver.saveParameter<Velocity<>, Lattice::NDIM>(timestep);
        }

        // Evolve by one timestep
        lbm.evolve();
    }
    return 0;
}
