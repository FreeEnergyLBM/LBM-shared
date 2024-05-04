#include <lbm.hh>

// This script simulates a Poiseuille flow in a channel driven by Zou-He boundary conditions

const int lx = 50;  // Channel length (+1)
const int ly = 20;  // Channel width (+4)

const int timesteps = 10000;
const int saveInterval = 10000;

const double pDiff = 1e-4;

using Lattice = LatticeProperties<ParallelX<1>, lx, ly>;

int initBoundaries(const int k) {
    int y = computeY<Lattice>(k);
    int x = computeX<Lattice>(k);
    if (y < 2 || y >= ly - 2) return 1;
    if (x == 0 || x == lx - 1) return 2;
    return 0;
}

double initDensity(int k) {
    int x = computeX<Lattice>(k);
    double t = x / (lx - 1.0);
    double rhoIn = 1 + 3 * pDiff / 2;
    double rhoOut = 1 - 3 * pDiff / 2;
    return rhoIn * (1 - t) + rhoOut * t;
}

int main(int argc, char **argv) {
    mpi.init();

    // Define the model
    using Trait = DefaultTraitFlowField<Lattice>::AddBoundary<ZouHeDensity>;
    FlowField<Lattice, Trait> model;

    // Set up geometry
    Geometry<Lattice>::initialiseBoundaries(initBoundaries);
    model.template getBoundary<BounceBack>().setNodeID({1});
    model.template getBoundary<ZouHeDensity>().setNodeID({2});
    model.setCollideID({0, 2});

    // Initialise fluid
    Density<>::template set<Lattice>(initDensity);

    SaveHandler<Lattice> saver("data/");
    saver.maskSolid();

    // Main loop
    Algorithm lbm(model);
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        if (timestep % saveInterval == 0) {
            print("Saving at timestep:", timestep);
            saver.saveVTK(timestep, Density<>::template getInstance<Lattice>(),
                          Velocity<>::template getInstance<Lattice, Lattice::NDIM>());
        }
        lbm.evolve();
    }

    return 0;
}
