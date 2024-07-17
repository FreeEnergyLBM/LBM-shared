#include <math.h>

#include <lbm.hh>

// This script simulates two phase flow through a porous structure in 3D.
// This example shows how to:
//      (i) pass the geometry to the lattice,
//      (ii) use saveBoundariesVTK to save the geometry, and
//      (iii) use saveVTK to save the simulation data.

const int lx = 250;  // Size of domain in x direction
const int ly = 250;  // Size of domain in y direction
const int lz = 10;   // Size of domain in z direction

const int timesteps = 5000000;   // Number of iterations to perform
const int saveInterval = 10000;  // Interval to save global data

// Parameters to control the surface tension and width of the diffuse interface
// Surface tension (in lattice units) = sqrt(8*kappa*A/9)
// Interface width (in lattice units) = sqrt(kappa/A)
// For reference, the model used is in section 9.2.2 of The Lattice Boltzmann Method: Principles and Practice, T. Kruger
// et al. (2016)
const double binaryA = 4.17e-3;
const double binaryKappa = 2 * binaryA;
const double contactAngle = 45;  // Contact angle of the liquid on the solid

// Relaxation times of each component. tau1 corresponds to phi=1.0, tau2 corresponds to phi=-1.0
// Viscosity (in lattice units) = 1.0/3.0 * (tau - 0.5)
double tau1 = 0.8;
double tau2 = 0.6;

const double inletOrderParameter = -1.0;  // Order parameter value at the inlet

// Set up the lattice, including the resolution and data/parallelisation method
using Lattice = LatticeProperties<ParallelX<1>, lx, ly, lz>;

void unzipGeometry() {
    if (mpi.rank == 0) {
        int result = system("unzip -o geometry.zip");
        if (result != 0) throw std::runtime_error("Failed to unzip geometry");
    }
    mpi.barrier();
}

// Function used to define the geometry
std::vector<int> initSolid() {
    // Read in the geometry file
    unzipGeometry();
    auto boundaries = loadTxt<int>("geometry.dat");
    // Define any additional boundaries
    for (int x = 0; x < lx; x++) {
        for (int y = 0; y < ly; y++) {
            for (int z = 0; z < lz; z++) {
                int k = x * ly * lz + y * lz + z;
                if (boundaries[k] != 0)
                    continue;
                else if (y <= 1 || y >= ly - 2 || x <= 1 || x >= lx - 2)
                    boundaries[k] = 1;
                else if (x == 2)
                    boundaries[k] = 3;
                else if (x == lx - 3)
                    boundaries[k] = 4;
            }
        }
    }
    return boundaries;
}

// Function used to define the fluid
// Here we set a tanh transition in the y direction. This should match the equilibrium profile.
double initFluid(int k) {
    int x = computeXGlobal<Lattice>(k);
    if (x == 2) return inletOrderParameter;
    return tanh((x - lx / 4.) / (sqrt(2 * binaryKappa / binaryA)));
}

using ConvectOP = ConvectParameterBoundary<OrderParameter<>, OrderParameterOld<>>;

int main(int argc, char **argv) {
    mpi.init();

    using TraitFlowFieldBinary = DefaultTraitFlowFieldBinary<Lattice>::SetCollisionOperator<MRT>::SetBoundary<
        std::tuple<VelocityInflow, Convective>, std::tuple<BounceBack>>;
    using TraitBinary = DefaultTraitBinary<Lattice>::template SetProcessor<
        std::tuple<ConvectOP>, std::tuple<GradientsMultiStencil<OrderParameter<>, CentralXYZ, LaplacianCentral>>,
        std::tuple<ChemicalPotentialCalculatorBinary, CubicWetting,
                   SetParameterOld<OrderParameter<>, OrderParameterOld<>>>>::
        SetBoundary<std::tuple<DirichletOrderParameter, Convective>, std::tuple<BounceBack>>;

    // Define the models to be used
    FlowFieldBinary<Lattice, TraitFlowFieldBinary>
        flowFieldModel;  // Flowfield (navier stokes solver) that can be used with the binary model
    Binary<Lattice, TraitBinary> componentSeparationModel;  // Binary model with hybrid equilibrium and forcing term

    // Pass the relaxation times to each model
    flowFieldModel.setTau1(tau1);
    flowFieldModel.setTau2(tau2);

    componentSeparationModel.setTau1(tau1);
    componentSeparationModel.setTau2(tau2);
    componentSeparationModel.setA(binaryA);

    std::cout << "dynamic viscosity ratio: " << 1.0 / 3.0 * (tau2 - 0.5) / (1.0 / 3.0 * (tau1 - 0.5)) << "\n";

    // Pass the surface tension/interface width parameters to the relevant preprocessor
    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setA(binaryA);
    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setKappa(binaryKappa);

    // Define the solid boundaries
    Geometry<Lattice>::initialiseBoundaries(initSolid());
    flowFieldModel.getBoundary<BounceBack>().setNodeID(1);
    componentSeparationModel.getBoundary<BounceBack>().setNodeID(1);

    std::vector<double> u(3);
    u[0] = -0.001;
    u[1] = 0.;
    u[2] = 0.;
    flowFieldModel.getBoundary<VelocityInflow>().setNodeID(3);
    flowFieldModel.getBoundary<VelocityInflow>().setWallVelocity(u);

    flowFieldModel.getBoundary<Convective>().setNodeID(4);

    componentSeparationModel.getBoundary<DirichletOrderParameter>().setNodeID(3);
    componentSeparationModel.getBoundary<DirichletOrderParameter>().setInterfaceVal(inletOrderParameter);

    componentSeparationModel.getBoundary<Convective>().setNodeID(4);
    componentSeparationModel.getProcessor<ConvectOP>().setNodeID(4);

    componentSeparationModel.getProcessor<CubicWetting>().setNodeID(1);
    componentSeparationModel.getProcessor<CubicWetting>().setThetaDegrees(contactAngle);
    componentSeparationModel.getProcessor<CubicWetting>().setAlpha(sqrt(2));
    componentSeparationModel.getProcessor<CubicWetting>().setNeutralWetLayerThickness(10);
    componentSeparationModel.getProcessor<CubicWetting>().useSinglePhaseCheck = true;

    // Initialise the fluid using the function above
    OrderParameter<>::set<Lattice>(initFluid);

    // Algorithm creates an object that can run our chosen LBM model
    Algorithm lbm(flowFieldModel, componentSeparationModel);

    // Set up the handler object for saving data
    SaveHandler<Lattice> saver("data/");
    saver.maskSolid();

    saver.saveBoundariesVTK(0);

    // Perform the main LBM loop
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        // Save the desired parameters, producing a binary file for each.
        if (timestep % saveInterval == 0) {
            saver.saveVTK(timestep, Density<>::template getInstance<Lattice>(),
                          OrderParameter<>::template getInstance<Lattice>(),
                          Velocity<>::template getInstance<Lattice, Lattice::NDIM>());

            std::cout << "Saving at timestep " << timestep << "." << std::endl;
        }

        // Evolve by one timestep
        lbm.evolve();
    }

    std::cout << "Simulation complete." << std::endl;
    return 0;
}
