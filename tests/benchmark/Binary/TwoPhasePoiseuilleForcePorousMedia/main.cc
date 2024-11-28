#include <lbm.hh>

// This script simulates a Poiseuille flow in a porous channel driven by an external force.
// The benchmark is based on the analytical solution presented in 10.1103/PhysRevE.66.036304.

const int lx = 100;  // Size of domain in x direction
const int ly = 5;    // Size of domain in y direction

const int timesteps = 5000;    // Number of iterations to perform
const int saveInterval = 1000;  // Interval to save global data

const double porosity1 = 0.2;         // Porosity of the domain
const double porosity2 = 0.7;         // Porosity of the domain
const double permeability1 = 7.5;  // Permeability of the domain
const double permeability2 = 2286.66;  // Permeability of the domain

// Irreducible wetting and non-wetting saturations
double Sw_ir = 0.0;
double Snw_ir = 0.0;

// Capillary pressure parameters
double P_0 = 1.0E-2;
double P_c_inf = 10 * P_0;
double lambda = 2.0;

// Relaxation times of each component.
double tau1 = 1.0;

using Lattice = LatticeProperties<ParallelX<1>, lx, ly>;

// Function used to initialise the solid (1) and fluid (0)
int initSolid(int k) {
    int y = computeY(ly, 1, k);
    if (y == 0 || y == ly - 1) {
        return 0;  // fluid (for now)
    } else {
        return 0;  // fluid
    }
}

// Function used to initialise the porosity
double initPorosity(const int k) {
    int x = computeXGlobal<Lattice>(k);
    if (x < lx / 2) {
        return porosity1;
    } else {
        return porosity2;
    }
}

// Function used to initialise the permeability
double initPermeability(const int k) {
    int x = computeXGlobal<Lattice>(k);
    if (x < lx / 2) {
        return permeability1;
    } else {
        return permeability2;
    }
}

// Function used to initialise the liquid (1) and gas (-1)
double initSaturation(int k) {
    int x = computeXGlobal<Lattice>(k);  // global function used because the x direction is split among the processors

    if (x < lx / 2) {
        return 0.2;
    } else {
        return 0.2;
    }
}

int main(int argc, char **argv) {
    mpi.init();

    // Set up the model
    using WettingPorousPoiseuilleTrait =
        DefaultTraitMultiphasePorousFlowField<Lattice>::SetCollisionOperator<SRT>::SetProcessor<
            std::tuple<UpdateRelPermAndPc>, std::tuple<GradientsMultiStencil<CapillaryPressure<>, CentralXYZ>>>::
            AddForce<CapillaryForcePorous<GuoMultiphasePorous, Gradient, WettingRelativePermeability<>, Velocity<>>>;

    MultiphasePorousFlowField<Lattice, WettingPorousPoiseuilleTrait> model;

    // Set the relaxation time
    model.setTau(tau1);

    // set irriducible wetting saturation for the model
    model.setSwIr(Sw_ir);

    // Set the solid
    Geometry<Lattice>::initialiseBoundaries(initSolid);
    model.template getBoundary<BounceBack>().setNodeID(1);

    // Set the porosity
    Porosity<>::set<Lattice>(initPorosity);

    // Set the permeability
    Permeability<>::set<Lattice>(initPermeability);

    // Set the saturation
    Saturation<>::set<Lattice>(initSaturation);

    // Set the irreducible wetting and non-wetting saturations
    model.template getProcessor<UpdateRelPermAndPc>().setSwIr(Sw_ir);
    model.template getProcessor<UpdateRelPermAndPc>().setSnwIr(Snw_ir);

    // Set the capillary pressure parameters
    model.template getProcessor<UpdateRelPermAndPc>().setP0(P_0);
    model.template getProcessor<UpdateRelPermAndPc>().setPcInf(P_c_inf);
    model.template getProcessor<UpdateRelPermAndPc>().setLambda(lambda);

    // Set the Leverett J-function parameters
    model.template getProcessor<UpdateRelPermAndPc>().setK0(permeability1);
    model.template getProcessor<UpdateRelPermAndPc>().setPhi0(porosity1);

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
            saver.saveVTK(timestep, Saturation<>::template getInstance<Lattice>(),
                          CapillaryPressure<>::template getInstance<Lattice>(),
                          WettingRelativePermeability<>::template getInstance<Lattice>(),
                          Velocity<>::template getInstance<Lattice, Lattice::NDIM>());
        }

        // Evolve by one timestep
        lbm.evolve();
    }

    std::cout << "Simulation complete." << std::endl;
    return 0;
}