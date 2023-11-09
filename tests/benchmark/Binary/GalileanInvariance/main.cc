#include <lbm.hh>

// This script simulates a droplet with bulk fluid movement

const int timesteps = 10000; // Number of iterations to perform
const int saveInterval = 10000; // Interval to save global data

const int lx = 60; // Size of domain in x direction
const int ly = 60; // Size of domain in y direction
const int radius = 20; // Droplet radius


using Lattice = LatticeProperties<NoParallel, lx, ly>;


double initFluid(const int k) {
    int x = computeXGlobal<Lattice>(k);
    int y = computeY(ly, 1, k);
    int r2 = pow(x-lx/2, 2) + pow(y-ly/2, 2);
    if (r2 < radius*radius) return 1;
    else return -1;
}

//User defined function to define some fluid initialisation (optional)
bool velocityLocation(const int k) {
    return true;
}

int main(int argc, char **argv) {
    mpi.init();

    // Set up the model
    FlowFieldBinary<Lattice> model1;
    Binary<Lattice> model2;

    model2.getPreProcessor<ChemicalPotentialCalculatorBinary>().setA(0.001);
    model2.getPreProcessor<ChemicalPotentialCalculatorBinary>().setKappa(0.001);

    OrderParameter<>::set<Lattice>(initFluid);

    Velocity<>::set<Lattice,2,0>([](const int k){ return 0.001; });
    Velocity<>::set<Lattice,2,1>([](const int k){ return 0.001; });

    // Create save handler
    ParameterSave<Lattice> saver("data/");
    saver.SaveHeader(timesteps, saveInterval);

    // Initialise
    Algorithm lbm(model1, model2);

    // Main loop
    for (int timestep=0; timestep<=timesteps; timestep++) {
        if (timestep%saveInterval==0) {
          saver.SaveParameter<Density<>>(timestep);
          saver.SaveParameter<OrderParameter<>>(timestep);
          saver.SaveParameter<Velocity<>,Lattice::NDIM>(timestep);
        }
        lbm.evolve();
    }

    return 0;
}
