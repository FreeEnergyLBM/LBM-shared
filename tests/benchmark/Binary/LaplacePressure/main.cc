#include <lbm.hh>

// This script simulates a stationary droplet

const int timesteps = 10000;     // Number of iterations to perform
const int saveInterval = 10000;  // Interval to save global data

const int lx = 60;      // Size of domain in x direction
const int ly = 60;      // Size of domain in y direction
const int radius = 20;  // Droplet radius

using Lattice = LatticeProperties<NoParallel, lx, ly>;

double initFluid(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    const int rr2 = pow(xx - lx / 2.0, 2) + pow(yy - ly / 2.0, 2);
    if (rr2 < radius * radius)
        return 1;
    else
        return -1;
}

int main(int argc, char **argv) {
    mpi.init();

    // Set up the model
    FlowFieldBinary<Lattice> model1;
    Binary<Lattice> model2;

    model2.template getProcessor<ChemicalPotentialCalculatorBinary>().setA(0.00015);
    model2.template getProcessor<ChemicalPotentialCalculatorBinary>().setKappa(0.0003);

    // Initialise fluid
    OrderParameter<>::template set<Lattice>(initFluid);

    // Create save handler
    SaveHandler<Lattice> saver("data/");
    saver.saveHeader(timesteps, saveInterval);

    // Main loop
    Algorithm lbm(model1, model2);
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        if (timestep % saveInterval == 0) {
            saver.template saveParameter<Density<>>(timestep);
            saver.template saveParameter<OrderParameter<>>(timestep);
        }
        lbm.evolve();
    }

    return 0;
}
