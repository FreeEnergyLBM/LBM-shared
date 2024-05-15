#include <lbm.hh>

// This script simulates 3D capillary rise.

const int lx = 40;   // Size of domain in x direction
const int ly = 100;  // Size of domain in y direction
const int lz = 40;   // Size of domain in z direction

const double capillaryRadius = lx / 8;  // Radius of the capillary

std::vector<double> gravity = {0.0, -1.0e-5};  // Gravitational force

const int timesteps = 10000000;  // Number of iterations to perform
const int saveInterval = 10000;  // Interval to save global data

// Surface tension (in lattice units) = sqrt(8 * kappa * A / 9)
// Interface width (in lattice units) = sqrt(kappa / A)
const double binaryA = 1.0e-3;
const double binaryKappa = 2 * binaryA;
const double contactAngle = 30;  // Contact angle of the liquid (phi = 1) on the solid

// Relaxation times of each component. tau1 corresponds to phi = 1.0, tau2 corresponds to phi = -1.0
double tau1 = 0.621;
double tau2 = 0.509;

// Set up the lattice, including the resolution and data/parallelisation method
using Lattice = LatticeProperties<ParallelX<1>, lx, ly, lz>;

// Function used to define the solid geometry
int initSolid(const int k) {
    int x = computeXGlobal<Lattice>(k);
    int y = computeY(ly, lz, k);
    int z = computeZ(ly, lz, k);

    double distxz = sqrt(pow(x - lx / 2, 2) + pow(z - lz / 2, 2));

    if (y <= 1 || y >= ly - 2)
        return 1;
    else if (y >= ly * 0.05 && y <= ly * 0.95) {
        if (distxz >= capillaryRadius && distxz < capillaryRadius + 2)
            return 1;
        else
            return 0;
    }
    return 0;
}

// Function used to initialise the liquid (1) and gas (-1)
double initFluid(int k) {
    int y = computeY(ly, lz, k);
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
    componentSeparationModel.setA(binaryA);

    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setA(binaryA);
    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setKappa(binaryKappa);

    // Set the boundaries
    Geometry<Lattice>::initialiseBoundaries(initSolid);

    flowFieldModel.getBoundary<BounceBack>().setNodeID(1);
    componentSeparationModel.getBoundary<BounceBack>().setNodeID(1);

    componentSeparationModel.getProcessor<CubicWetting>().setNodeID(1);
    componentSeparationModel.getProcessor<CubicWetting>().setThetaDegrees(contactAngle);
    componentSeparationModel.getProcessor<CubicWetting>().setAlpha(sqrt(2));

    // Set the force acting on the fluid: first argument is the force vector, second argument is the component to which
    // it applies. 0 corresponds to the liquid component (phi = 1), 1 corresponds to the gas component (phi = -1).
    flowFieldModel.getForce<BodyForce<>>().setForce(gravity, 0);

    // Initialise the liquid and gas
    OrderParameter<>::set<Lattice>(initFluid);

    // Algorithm creates an object that combines the lattice models and runs them in order
    Algorithm lbm(flowFieldModel, componentSeparationModel);

    // Set up the handler object for saving data
    SaveHandler<Lattice> saver("data/");
    saver.maskSolid();

    saver.saveBoundariesVTK(0);
    saver.saveHeader(timesteps, saveInterval);  // Create a header with lattice information (lx, ly, lz, NDIM (2D or
                                                // 3D), timesteps, saveInterval)

    // Perform the main LBM loop
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        if (timestep % saveInterval == 0) {
            saver.saveVTK(timestep, Density<>::template getInstance<Lattice>(),
                          OrderParameter<>::template getInstance<Lattice>(),
                          Velocity<>::template getInstance<Lattice, Lattice::NDIM>());

            saver.saveParameter<OrderParameter<>>(timestep);
            saver.saveParameter<Velocity<>, Lattice::NDIM>(timestep);
            std::cout << "Saving at timestep " << timestep << "." << std::endl;
        }
        // Evolve by one timestep
        lbm.evolve();
    }

    return 0;
}
