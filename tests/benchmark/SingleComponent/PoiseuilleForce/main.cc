#include <lbm.hh>

// This script simulates a Poiseuille flow in a channel driven by an external force

const int lx = 2;   // Size of domain in x direction
const int ly = 50;  // Size of domain in y direction

const int timesteps = 10000;     // Number of iterations to perform
const int saveInterval = 10000;  // Interval to save global data

const double force = 1e-6;  // Driving force, equivalent to the pressure gradient

using Lattice = LatticeProperties<NoParallel, lx, ly>;

// Set a solid at the top and bottom
int initSolid(const int k) {
    int y = computeY(ly, 1, k);
    if (y <= 1 || y >= ly - 2)
        return 1;
    else
        return 0;
}

int main(int argc, char **argv) {
    mpi.init();

    // Set up the model
    using PoiseuilleTrait = DefaultTraitFlowField<Lattice>::AddForce<BodyForce<>>;
    FlowField<Lattice, PoiseuilleTrait> model;
    model.template getForce<BodyForce<>>().setForce({force, 0, 0});

    // Set the solid
    Geometry<Lattice>::initialiseBoundaries(initSolid);
    model.template getBoundary<BounceBack>().setNodeID(1);

    // Create save handler
    SaveHandler<Lattice> saver("data/");
    saver.saveHeader(timesteps, saveInterval);

    // Initialise
    Algorithm lbm(model);

    // Main loop
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        if (timestep % saveInterval == 0) saver.template saveParameter<Velocity<>, Lattice::NDIM>(timestep);
        lbm.evolve();
    }

    return 0;
}
