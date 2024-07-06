#include <lbm.hh>

#ifdef MPIPARALLEL
#define PARALLEL
#endif
#include <minim.h>

const int lx = 100;
const int ly = 100;

const int timesteps = 5000;
const int saveInterval = 1000;

const double channelWidth = 40;
const double force = 1e-6;

std::vector<bool> initialiseSolid() {
    const double PI = acos(-1);
    std::vector<bool> solid(lx * ly);
    for (int x = 0; x < lx; x++) {
        for (int y = 0; y < ly; y++) {
            double phase = 2 * PI * x / (lx - 1);
            double y1 = 0.5 * (cos(phase) + 1) * (ly - 1 - channelWidth);
            double y2 = y1 + channelWidth;
            if (y <= y1 || y > y2) solid[x * ly + y] = true;
        }
    }
    return solid;
}

std::vector<double> initialiseFluid(std::vector<double> solid) {
    std::vector<double> data(3 * lx * ly);
    for (int x = 0; x < lx; x++) {
        for (int y = 0; y < ly; y++) {
            int i = x * ly + y;

            // Initialise the liquid
            if (x < 0.5 * (lx - 1)) {
                data[3 * i + 1] = 1;  // Liquid
            } else {
                data[3 * i + 2] = 1;  // Gas
            }

            // Account for the solid phase
            data[3 * i + 0] = solid[3 * i];
            data[3 * i + 1] *= (1 - solid[3 * i]);
            data[3 * i + 2] *= (1 - solid[3 * i]);
        }
    }
    return data;
}

int main(int argc, char** argv) {
#ifdef MPIPARALLEL
    mpi.init();
    minim::mpi.init(MPI_COMM_WORLD);
#endif

    // ENERGY MINIMISATION

    // Set up the potential
    minim::PhaseField potential;
    potential.setNFluid(3);
    potential.setGridSize({lx, ly, 1});
    potential.setDensityConstraint("hard");
    potential.setSolid([](int x, int y, int z) { return (x == 0) || (x == lx - 1); });

    // Diffuse-solid initialisation
    print("Initialising solid");
    std::vector<double> solid = potential.diffuseSolid(initialiseSolid());
    potential.setFixFluid(0);
    potential.setVolumeFixed(true);

    // Minimise the liquid
    print("Minimising liquid");
    auto state = minim::State(potential, initialiseFluid(solid));
    auto logFn = [](int iter, minim::State& state) {
        if (iter % 100 == 0) print(iter, ": Energy =", state.energy());
    };
    auto minimiser = minim::Lbfgs().setMaxIter(1000);
    std::vector<double> minimum = minimiser.minimise(state, logFn);
    print("Converged after", minimiser.iter, "iterations");

    std::ofstream file("minimum.txt");
    for (int i = 0; i < lx * ly * 3; i++) {
        file << minimum[i] << std::endl;
    }
    file.close();

    // LATTICE BOLTZMANN

    using Lattice = LatticeProperties<ParallelX<1>, lx, ly>;

    // Initialise based on the minimised coordinates
    print("Initialising LBM");
    for (auto [x, y, z, k] : RangeXYZK<Lattice>()) {
        int kGlobal = x * ly + y;
        OrderParameter<0>::template initialise<Lattice>(k, minimum[3 * kGlobal]);
        OrderParameter<1>::template initialise<Lattice>(k, minimum[3 * kGlobal + 1]);
    }

    // Models
    FlowField<Lattice> flowFieldModel;
    TernaryLee<Lattice, 0> ternaryModel1;
    TernaryLee<Lattice, 1> ternaryModel2;
    Algorithm lbm(flowFieldModel, ternaryModel1, ternaryModel2);

    // LBM loop
    print("Running LBM");
    SaveHandler<Lattice> saver("data/");
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        if (timestep % saveInterval == 0) {
            print("Saving at timestep", timestep);
            saver.saveVTK(timestep, Velocity<>::template getInstance<Lattice, 2>(),
                          OrderParameter<0>::template getInstance<Lattice>(),
                          OrderParameter<1>::template getInstance<Lattice>());
        }

        lbm.evolve();
    }

    return 0;
}
