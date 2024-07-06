#include <lbm.hh>

// This script simulatee single phase solute transport in a 2D pipe.

const int lx = 300;  // Size of domain in x direction
const int ly = 100;  // Size of domain in y direction

const int timesteps = 500000;    // Number of iterations to perform
const int saveInterval = 10000;  // Interval to save global data

const double inletC = 1.0;         // concentration at the inlet
const double diffusivity = 0.001;  // diffusivity of the solute

// Set up the lattice, including the resolution and data/parallelisation method
using Lattice = LatticeProperties<ParallelX<1>, lx, ly>;

// Function used to define the solid geometry
// Here we set a solid at the top and bottom
int initSolid(const int k) {
    int x = computeXGlobal<Lattice>(k);
    int y = computeY(ly, 1, k);
    if (y <= 1 || y >= ly - 2) return 1;
    if (x <= 1 || x >= lx - 2)
        return 1;
    else if (x == 2)
        return 3;
    else if (x == lx - 3)
        return 4;
    return 0;
}

double initSolute(int k) {
    int x = computeXGlobal<Lattice>(k);
    if (x == 2) return inletC;
    return 0.0;
}

using ConvectSolute = ConvectParameterBoundary<Solute<>, SoluteOld<>>;

int main(int argc, char **argv) {
    mpi.init();

    // We use the 'FlowField' LBM model, which is the standard Navier-Stokes solver
    // We need to modify the traits of the model to include a body force.
    using TraitFlowField =
        DefaultTraitFlowField<Lattice>::SetCollisionOperator<MRT>::SetBoundary<std::tuple<BounceBack>,
                                                                               std::tuple<VelocityInflow, Convective>>;
    using ADTrait = DefaultTraitAdvectionDiffusion<Lattice>::template SetProcessor<
        std::tuple<ConvectSolute>, std::tuple<SetParameterOld<Solute<>, SoluteOld<>>>>::
        SetBoundary<std::tuple<BounceBack>, std::tuple<Dirichlet, Convective>>;

    FlowField<Lattice, TraitFlowField> flowFieldModel;
    AdvectionDiffusion<Solute<>, Lattice, ADTrait> ADModel;

    // Set the diffusivity of the solute
    ADModel.setDiffusivity(diffusivity);

    // Define the solid boundaries
    Geometry<Lattice>::initialiseBoundaries(initSolid);
    flowFieldModel.getBoundary<BounceBack>().setNodeID(1);
    ADModel.getBoundary<BounceBack>().setNodeID(1);

    std::vector<double> u(3);
    u[0] = -0.001;
    u[1] = 0.;
    u[2] = 0.;
    flowFieldModel.getBoundary<VelocityInflow>().setNodeID(3);
    flowFieldModel.getBoundary<VelocityInflow>().setWallVelocity(u);

    flowFieldModel.getBoundary<Convective>().setNodeID(4);

    ADModel.getBoundary<Dirichlet>().setNodeID(3);
    ADModel.getBoundary<Dirichlet>().setInterfaceVal(inletC);

    ADModel.getBoundary<Convective>().setNodeID(4);
    ADModel.getProcessor<ConvectSolute>().setNodeID(4);

    Solute<>::set<Lattice>(initSolute);

    // Set up the handler object for saving data
    SaveHandler<Lattice> saver("data/");
    saver.maskSolid();

    saver.saveBoundariesVTK(0);

    // Algorithm creates an object that can run our chosen LBM model
    Algorithm lbm(flowFieldModel, ADModel);

    // Perform the main LBM loop
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        // Save the desired parameters, producing a binary file for each.
        if (timestep % saveInterval == 0) {
            std::cout << "Saving at timestep " << timestep << "." << std::endl;
            saver.saveVTK(timestep, Density<>::template getInstance<Lattice>(),
                          Solute<>::template getInstance<Lattice>(),
                          Velocity<>::template getInstance<Lattice, Lattice::NDIM>());
        }

        // Evolve by one timestep
        lbm.evolve();
    }

    return 0;
}
