#include <lbm.hh>

// This script simulates 3D Young-Laplace test. Check the Verification.xlsx for the expected results.

const int lx = 100;                  // Size of domain in x direction
const int ly = 100;                  // Size of domain in y direction
const int lz = 100;                  // Size of domain in z direction
const double dropRadius = 0.2 * lx;  // Radius to initialise the droplet

const int timesteps = 200000;    // Number of iterations to perform
const int saveInterval = 10000;  // Interval to save global data

// Relaxation times of each component. tau1 corresponds to phi = 1.0, and tau2 corresponds to phi = -1.0.
double tau1 = 0.6;
double tau2 = 0.6;

// Parameters to control the surface tension and width of the diffuse interface
// Surface tension (in lattice units) = sqrt(8 * kappa * A / 9)
// Interface width (in lattice units) = sqrt(kappa / A)
double A = 0.0015;
double kappa = A;

// Set up the lattice, including the resolution and data/parallelisation method
using Lattice = LatticeProperties<ParallelX<1>, lx, ly, lz>;

// Function used to initialise the liquid (1) and gas (-1)
double initFluid(int k) {
    int x = computeXGlobal<Lattice>(k);  // global function used because the x direction is split among the processors
    int y = computeY(ly, lz, k);
    int z = computeZ(ly, lz, k);

    double x0 = lx / 2.0;
    double y0 = ly / 2.0;
    double z0 = lz / 2.0;

    double distance = sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2));

    if (distance < dropRadius) {
        return 1;
    } else {
        return -1;
    }
}

// Modify the traits of the binary model to use MRT
using TraitFlowFieldBinary = DefaultTraitFlowFieldBinary<Lattice>::SetCollisionOperator<MRT>;

int main(int argc, char **argv) {
    mpi.init();

    // Define the models to be used
    FlowFieldBinary<Lattice, TraitFlowFieldBinary>
        flowFieldModel;  // Flowfield (navier stokes solver) that can be used with the binary model
    Binary<Lattice> componentSeparationModel;  // Binary model with hybrid equilibrium and forcing term

    // Pass the relaxation times to each model
    flowFieldModel.setTau1(tau1);
    flowFieldModel.setTau2(tau2);

    componentSeparationModel.setTau1(tau1);
    componentSeparationModel.setTau2(tau2);

    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setA(A);
    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setKappa(kappa);

    // Initialise the liquid and gas
    OrderParameter<>::set<Lattice>(initFluid);

    // Algorithm creates an object that combines the lattice models and runs them in order
    Algorithm lbm(flowFieldModel, componentSeparationModel);

    // Set up the handler object for saving data
    SaveHandler<Lattice> saver("data/");

    // Perform the main LBM loop
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        if (timestep % saveInterval == 0) {
            std::cout << "Saving at timestep " << timestep << "." << std::endl;
            saver.saveVTK(timestep, Density<>::template getInstance<Lattice>(),
                          OrderParameter<>::template getInstance<Lattice>());
        }
        lbm.evolve();
    }

    return 0;
}
