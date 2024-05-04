#include <lbm.hh>

// This script simulates capillary rise.

const int lx = 200;  // Size of domain in x direction
const int ly = 200;  // Size of domain in y direction

const double capillaryRadius = 0.05 * lx;

const double force = -2.5e-6;  // Driving force, equivalent to the pressure gradient

const int timesteps = 10000000;  // Number of iterations to perform
const int saveInterval = 10000;  // Interval to save global data

const double binaryA = 1.0e-3;
const double binaryKappa = 2 * binaryA;
const double contactAngle = 30;  // Contact angle of the liquid on the solid

// Relaxation times of each component. tau1 corresponds to phi=1.0, tau2 corresponds to phi=-1.0
// Viscosity (in lattice units) = 1.0/3.0 * (tau - 0.5)
double tau1 = 0.621;
double tau2 = 0.509;

// Set up the lattice, including the resolution and data/parallelisation method
using Lattice = LatticeProperties<ParallelX<1>, lx, ly>;

// Function used to define the solid geometry
int initSolid(const int k) {
    int x = computeXGlobal<Lattice>(k);
    int y = computeY(ly, 1, k);
    if (y <= 1 || y >= ly - 2)
        return 1;
    else if (y >= ly * 0.05 && y <= ly * 0.95) {
        if (abs(x - lx / 2) >= capillaryRadius && abs(x - lx / 2) <= capillaryRadius + 2)
            return 1;
        else
            return 0;
    }
    return 0;
}

// Function used to initialise the liquid (1) and gas (-1)
double initFluid(int k) {
    int y = computeY(ly, 1, k);
    return -tanh((y - ly / 2) / (sqrt(2 * binaryKappa / binaryA)));
}

// Modify the traits of the binary model to use MRT
using TraitFlowFieldBinary = DefaultTraitFlowFieldBinary<Lattice>::SetCollisionOperator<MRT>::AddForce<BodyForce<>>;

int main(int argc, char **argv) {
    mpi.init();

    // Define the models to be used
    FlowFieldBinary<Lattice, TraitFlowFieldBinary>
        flowFieldModel;  // Flowfield (navier stokes solver) that can be used with the binary model
    Binary<Lattice> componentSeparationModel;  // Binary model with hybrid equilibrium and forcing term

    // Set the relaxation times for the lattice models
    flowFieldModel.setTau1(tau1);
    flowFieldModel.setTau2(tau2);

    componentSeparationModel.setTau1(tau1);
    componentSeparationModel.setTau2(tau2);

    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setA(binaryA);
    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setKappa(binaryKappa);

    // Set the boundaries
    Geometry<Lattice>::initialiseBoundaries(initSolid);

    flowFieldModel.getBoundary<BounceBack>().setNodeID(1);
    componentSeparationModel.getBoundary<BounceBack>().setNodeID(1);

    componentSeparationModel.getProcessor<CubicWetting>().setNodeID(1);
    componentSeparationModel.getProcessor<CubicWetting>().setThetaDegrees(contactAngle);
    componentSeparationModel.getProcessor<CubicWetting>().setAlpha(sqrt(2));

    // Define the magnitude of the body force
    flowFieldModel.getForce<BodyForce<>>().setMagnitudeY(force);
    flowFieldModel.getForce<BodyForce<>>().activateGravityY();

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
            saver.saveVTK(timestep, Density<>::template getInstance<Lattice>(),
                          OrderParameter<>::template getInstance<Lattice>(),
                          Velocity<>::template getInstance<Lattice, Lattice::NDIM>());
        }
        lbm.evolve();
    }

    return 0;
}
