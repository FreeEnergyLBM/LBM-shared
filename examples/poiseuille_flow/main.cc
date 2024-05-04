#include <lbm.hh>

// This script simulates a Poiseuille flow in a channel driven by an external force

const int lx = 100;  // Size of domain in x direction
const int ly = 50;   // Size of domain in y direction

const int timesteps = 10000;    // Number of iterations to perform
const int saveInterval = 1000;  // Interval to save global data

const double force = 1e-6;  // Driving force, equivalent to the pressure gradient

// Function used to define the solid geometry
// Here we set a solid at the top and bottom
int initSolid(const int k) {
    int y = computeY(ly, 1, k);
    if (y <= 1 || y >= ly - 2)
        return 1;
    else
        return 0;
}

int main(int argc, char **argv) {
    mpi.init();

    // Set up the lattice, including the resolution and data/parallelisation method
    using Lattice = LatticeProperties<ParallelX<1>, lx, ly>;

    // We use the 'FlowField' LBM model, which is the standard Navier-Stokes solver
    // We need to modify the traits of the model to include a body force.
    using PoiseuilleTrait = DefaultTraitFlowField<Lattice>::AddForce<BodyForce<>>;
    FlowField<Lattice, PoiseuilleTrait> model;

    // Define the magnitude of the body force
    model.getForce<BodyForce<>>().setMagnitudeX(force);

    // Define the solid boundaries
    Geometry<Lattice>::initialiseBoundaries(initSolid);
    model.getBoundary<BounceBack>().setNodeID(1);

    // Set up the handler object for saving data
    SaveHandler<Lattice> saver("data/");
    saver.saveHeader(timesteps, saveInterval);  // Create a header with lattice information (lx, ly, lz, NDIM (2D or
                                                // 3D), timesteps, saveInterval)

    // Algorithm creates an object that can run our chosen LBM model
    Algorithm lbm(model);

    // Perform the main LBM loop
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        if (timestep % saveInterval == 0) {
            // Use the save handler to save the solid and the velocity
            saver.saveParameter<Velocity<>, Lattice::NDIM>(timestep);
            std::cout << "Saving at timestep " << timestep << "." << std::endl;
        }
        lbm.evolve();
    }

    return 0;
}
