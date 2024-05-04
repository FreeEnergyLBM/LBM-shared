#include <lbm.hh>

// This script simulates a flow in a channel with two liquids driven by an external force

const int lx = 2;    // Size of domain in x direction
const int ly = 100;  // Size of domain in y direction

const int timesteps = 50000;     // Number of iterations to perform
const int saveInterval = 10000;  // Interval to save global data

const double force = 1e-6;  // Driving force, equivalent to the pressure gradient

using Lattice = LatticeProperties<NoParallel, lx, ly>;

int initSolid(const int k) {
    int y = computeY(ly, 1, k);
    if (y <= 1 || y >= ly - 2)
        return 1;
    else
        return 0;
}

double initFluid(const int k) {
    int yy = computeY(ly, 1, k);
    if (yy < 0.5 * (ly - 1))
        return 1;
    else
        return -1;
}

int main(int argc, char **argv) {
    mpi.init();

    // Set up the model
    using PoiseuilleTrait = DefaultTraitFlowField<Lattice>::AddForce<BodyForce<>>;
    FlowFieldBinary<Lattice, PoiseuilleTrait> model1;
    Binary<Lattice> model2;

    model1.template getForce<BodyForce<>>().setForce({force, 0, 0});
    model2.setTau1(0.55);
    model2.setTau2(1.0);

    // Set the solid
    Geometry<Lattice>::initialiseBoundaries(initSolid);
    model1.template getBoundary<BounceBack>().setNodeID(1);
    model2.template getBoundary<BounceBack>().setNodeID(1);

    // Initialise the liquid
    OrderParameter<>::template set<Lattice>(initFluid);

    // Create save handler
    SaveHandler<Lattice> saver("data/");
    saver.saveHeader(timesteps, saveInterval);

    // Initialise
    Algorithm lbm(model1, model2);

    // Main loop
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        if (timestep % saveInterval == 0) saver.template saveParameter<Velocity<>, Lattice::NDIM>(timestep);
        lbm.evolve();
    }

    return 0;
}
