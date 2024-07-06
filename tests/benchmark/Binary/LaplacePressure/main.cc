#include <lbm.hh>

// This script simulates a stationary droplet

const int timesteps = 5000;     // Number of iterations to perform
const int saveInterval = 1000;  // Interval to save global data

const int lx = 60;      // Size of domain in x direction
const int ly = 60;      // Size of domain in y direction
const int radius = 20;  // Droplet radius

// Phase-field parameters
const double pfA = 0.00015;
const double pfKappa = 0.0003;

using Lattice = LatticeProperties<NoParallel, lx, ly>;

double initOrderParam(const int k) {
    double widthSqrt2 = sqrt(2 * pfKappa / pfA);
    auto [x, y, z] = computeXYZ<Lattice>(k);
    double r = sqrt(pow(x - lx / 2.0, 2) + pow(y - ly / 2.0, 2));
    return tanh((radius - r) / widthSqrt2);
}

int main(int argc, char **argv) {
    mpi.init();

    // Set up the model
    FlowFieldBinary<Lattice> model1;
    Binary<Lattice> model2;

    model2.template getProcessor<ChemicalPotentialCalculatorBinary>().setA(pfA);
    model2.template getProcessor<ChemicalPotentialCalculatorBinary>().setKappa(pfKappa);

    // Initialise fluid
    OrderParameter<>::template set<Lattice>(initOrderParam);

    // Create save handler
    SaveHandler<Lattice> saver("data/");
    saver.saveHeader(timesteps, saveInterval);

    // Main loop
    Algorithm lbm(model1, model2);
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        if (timestep % saveInterval == 0) {
            print("Saving at timestep:", timestep);
            saver.template saveParameter<Density<>>(timestep);
            saver.template saveParameter<OrderParameter<>>(timestep);
        }
        lbm.evolve();
    }

    return 0;
}
