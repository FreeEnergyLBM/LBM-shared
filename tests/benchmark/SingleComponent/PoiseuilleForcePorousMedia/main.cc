#include <lbm.hh>

// This script simulates a Poiseuille flow in a porous channel driven by an external force.
// The benchmark is based on the analytical solution presented in 10.1103/PhysRevE.66.036304.

const int lx = 84;  // Size of domain in x direction
const int ly = 84;  // Size of domain in y direction

const int timesteps = 500;     // Number of iterations to perform
const int saveInterval = 100;  // Interval to save global data

const double porosity = 0.1;          // Porosity of the medium
const double permeability = 6.40E-2;  // Permeability of the medium

const double force = 2.17e-6;  // Driving force, equivalent to the pressure gradient
const double tau = 0.6;       // Relaxation time

using Lattice = LatticeProperties<NoParallel, lx, ly>;

// Set a solid at the top and bottom
int initSolid(const int k) {
    int y = computeY(ly, 1, k);
    if (y <= 1 || y >= ly - 2)
        return 1;
    else
        return 0;
}

double initPorosity(const int k) { return porosity; }
double initPermeability(const int k) { return permeability; }

int main(int argc, char **argv) {
    mpi.init();

    std::cout << "Darcy Number = " << permeability / ((ly - 4) * (ly - 4)) << std::endl;
    double nu = 1.0 / 3.0 * (tau - 0.5);
    double u_0 = force * permeability * (1 - 1 / cosh(sqrt(porosity / permeability) * (ly - 4) / 2.0)) / nu;
    std::cout << "Reynolds = " << u_0 * (ly - 4) / nu << std::endl;
    std::cout << "Analytical u_max = " << u_0 << std::endl;

    // Set up the model
    using PorousPoiseuilleTrait = DefaultTraitPorousFlowField<Lattice>::AddForce<BodyForcePorous<>>;
    PorousFlowField<Lattice, PorousPoiseuilleTrait> model;
    model.template getForce<BodyForcePorous<>>().setForce({force, 0, 0});

    // Set the relaxation time
    model.setTau(tau);

    // Set the solid
    Geometry<Lattice>::initialiseBoundaries(initSolid);
    model.template getBoundary<BounceBack>().setNodeID(1);

    // Set the porosity
    Porosity<>::set<Lattice>(initPorosity);

    // Set the permeability
    Permeability<>::set<Lattice>(initPermeability);

    // Initialise the model
    Algorithm lbm(model);

    // Set up the handler object for saving data
    SaveHandler<Lattice> saver("data/");
    saver.maskSolid();

    saver.saveBoundariesVTK(0);

    // Perform the main LBM loop
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        // Save the desired parameters, producing a binary file for each.
        if (timestep % saveInterval == 0) {
            std::cout << "Saving at timestep " << timestep << "." << std::endl;
            saver.saveVTK(timestep, Density<>::template getInstance<Lattice>(),
                          Porosity<>::template getInstance<Lattice>(), Permeability<>::template getInstance<Lattice>(),
                          Velocity<>::template getInstance<Lattice, Lattice::NDIM>());
        }

        // Evolve by one timestep
        lbm.evolve();
    }

    std::cout << "Simulation complete." << std::endl;
    return 0;
}