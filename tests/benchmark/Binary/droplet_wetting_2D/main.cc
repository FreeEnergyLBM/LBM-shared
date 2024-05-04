#include <lbm.hh>

// This script simulates a droplet on a flat surface with a given contact angle.

const int lx = 50;  // Size of domain in x direction
const int ly = 50;  // Size of domain in y direction

const int timesteps = 1000;    // Number of iterations to perform
const int saveInterval = 100;  // Interval to save global data

const double contactAngle = 30;      // Contact angle of the liquid on the solid
const double dropRadius = 0.2 * lx;  // Radius to initialise the droplet

// Set up the lattice, including the resolution and data/parallelisation method
using Lattice = LatticeProperties<ParallelX<1>, lx, ly>;

// Function used to define the solid geometry
int initSolid(const int k) {
    int y = computeY(ly, 1, k);
    if (y <= 1 || y >= ly - 2)
        return 1;
    else
        return 0;
}

// Function used to initialise the liquid (1) and gas (-1)
double initFluid(int k) {
    int x = computeXGlobal<Lattice>(k);  // global function used because the x direction is split among the processors
    int y = computeY(ly, 1, k);

    // Check if within droplet radius
    double y0 = 1.5;  // - dropRadius * cos(contactAngle*M_PI/180);
    double r2 = pow(x - lx / 2.0, 2) + pow(y - y0, 2);
    if (r2 < pow(dropRadius, 2))
        return 1;
    else
        return -1;
}

// Modify the traits of the binary model to use MRT
using TraitFlowFieldBinary = DefaultTraitFlowFieldBinary<Lattice>::SetCollisionOperator<MRT>;

int main(int argc, char **argv) {
    mpi.init();

    // Define the models to be used
    FlowFieldBinary<Lattice, TraitFlowFieldBinary>
        flowFieldModel;  // Flowfield (navier stokes solver) that can be used with the binary model
    Binary<Lattice> componentSeparationModel;  // Binary model with hybrid equilibrium and forcing term

    componentSeparationModel.setTau2(0.51);
    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setA(0.015);
    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setKappa(0.03);

    // Set the boundaries
    Geometry<Lattice>::initialiseBoundaries(initSolid);

    flowFieldModel.getBoundary<BounceBack>().setNodeID(1);
    componentSeparationModel.getBoundary<BounceBack>().setNodeID(1);

    componentSeparationModel.getProcessor<CubicWetting>().setNodeID(1);
    componentSeparationModel.getProcessor<CubicWetting>().setThetaDegrees(contactAngle);
    componentSeparationModel.getProcessor<CubicWetting>().setAlpha(sqrt(2));

    // Initialise the liquid and gas
    OrderParameter<>::set<Lattice>(initFluid);

    // Algorithm creates an object that combines the lattice models and runs them in order
    Algorithm lbm(flowFieldModel, componentSeparationModel);

    // Set up the handler object for saving data
    SaveHandler<Lattice> saver("data/");
    saver.maskSolid();

    // Perform the main LBM loop
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        if (timestep % saveInterval == 0) {
            std::cout << "Saving at timestep " << timestep << "." << std::endl;
            saver.saveDAT(timestep, Density<>::template getInstance<Lattice>(),
                          OrderParameter<>::template getInstance<Lattice>(),
                          Velocity<>::template getInstance<Lattice, Lattice::NDIM>());
        }
        lbm.evolve();
    }

    return 0;
}
