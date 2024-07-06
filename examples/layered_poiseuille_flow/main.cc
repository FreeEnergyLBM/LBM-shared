#include <math.h>

#include <lbm.hh>

// This script simulates a Poiseuille flow in a channel with two immiscible layers of fluid driven by an external force.
// You can modify tau1 and tau2 to see how the velocity profile changes. tau1 and tau2 should be less than 10 or the
// accuracy degrades. They should never be less than 0.5, and values very close to 0.5 are prone to instability.

const int lx = 100;  // Size of domain in x direction
const int ly = 100;  // Size of domain in y direction

const int timesteps = 10000;    // Number of iterations to perform
const int saveInterval = 1000;  // Interval to save global data

const double force = 1e-6;  // Driving force, equivalent to the pressure gradient

// Parameters to control the surface tension and width of the diffuse interface
// Surface tension (in lattice units) = sqrt(8*kappa*A/9)
// Interface width (in lattice units) = sqrt(kappa/A)
// For reference, the model used is in section 9.2.2 of The Lattice Boltzmann Method: Principles and Practice, T. Kruger
// et al. (2016)
double A = 0.00015;
double kappa = 0.0003;

// Relaxation times of each component. tau1 corresponds to phi=1.0, tau2 corresponds to phi=-1.0
// Viscosity (in lattice units) = 1.0/3.0 * (tau - 0.5)
double tau1 = 1.0;
double tau2 = 0.55;

// Function used to define the solid geometry
// Here we set a solid at the top and bottom, in the conditions that return 1;
int initSolid(const int k) {
    int y = computeY(ly, 1, k);
    if (y <= 1 || y >= ly - 2)
        return 1;
    else
        return 0;
}

// Function used to define the fluid
// Here we set a tanh transition in the y direction. This should match the equilibrium profile.
double initFluid(int k) {
    int y = computeY(ly, 1, k);
    return tanh((y - ly / 2.) / (sqrt(2 * kappa / A)));
}

int main(int argc, char **argv) {
    mpi.init();

    // Set up the lattice, including the resolution and data/parallelisation method
    using Lattice = LatticeProperties<ParallelX<1>, lx, ly>;

    // We need to modify the traits of the model to include a body force as an 'AddOn'.
    // We modify the default traits for the 'FlowFieldBinary' model, adding a bodyforce and setting the collision model
    // to MRT, which improves accuracy at higher viscosity ratios
    using TraitFlowField = DefaultTraitFlowFieldBinary<Lattice>::AddForce<BodyForce<>>::SetCollisionOperator<MRT>;

    // Define the models to be used
    FlowFieldBinary<Lattice, TraitFlowField>
        flowFieldModel;  // Flowfield (navier stokes solver) that can be used with the binary model
    Binary<Lattice> componentSeparationModel;  // Binary model with hybrid equilibrium and forcing term

    // Pass the relaxation times to each model
    flowFieldModel.setTau1(tau1);
    flowFieldModel.setTau2(tau2);
    componentSeparationModel.setTau1(tau1);
    componentSeparationModel.setTau2(tau2);

    // Pass the surface tension/interface width parameters to the relevant preprocessor
    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setA(A);
    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setKappa(kappa);

    // Define the magnitude of the body force
    flowFieldModel.getForce<BodyForce<>>().setForce({force, 0, 0});

    // Define the solid boundaries
    Geometry<Lattice>::initialiseBoundaries(initSolid);
    flowFieldModel.getBoundary<BounceBack>().setNodeID(1);
    componentSeparationModel.getBoundary<BounceBack>().setNodeID(1);

    // Initialise the fluid using the function above
    OrderParameter<>::set<Lattice>(initFluid);

    // Algorithm creates an object that can run our chosen LBM model
    Algorithm lbm(flowFieldModel, componentSeparationModel);

    // Set up the handler object for saving data
    SaveHandler<Lattice> saver("data/");
    saver.saveHeader(timesteps, saveInterval);  // Create a header with lattice information (lx, ly, lz, NDIM (2D or
                                                // 3D), timesteps, saveInterval)

    // Perform the main LBM loop
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        // Save the desired parameters, producing a binary file for each.
        if (timestep % saveInterval == 0) {
            std::cout << "Saving at timestep " << timestep << "." << std::endl;
            saver.saveParameter<OrderParameter<>>(timestep);
            saver.saveParameter<Velocity<>, Lattice::NDIM>(timestep);
        }

        // Evolve by one timestep
        lbm.evolve();
    }

    return 0;
}
