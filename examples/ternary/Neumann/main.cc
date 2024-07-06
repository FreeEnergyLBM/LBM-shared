#include "main.hh"

#include <chrono>
#include <cstdlib>
#include <thread>

int main(int argc, char **argv) {
#ifdef MPIPARALLEL
    mpi.init();

    initMPIBoundary<Lattice>();
#endif

    initParams("input.txt");

    auto ternary1 = initTernary<0>();
    auto ternary2 = initTernary<1>();
    auto pressure = initPressure<>();

    SaveHandler<Lattice> saver(datadir);

    OrderParameter<1>::set<Lattice, 1>(initFluid2);
    OrderParameter<1>::smooth<Lattice, 1>(4);
    OrderParameter<0>::set<Lattice, 1>(initFluid1);
    OrderParameter<0>::smooth<Lattice, 1>(4);
    // OrderParameter<1>::set<Lattice,1>(0.0);
    ChemicalPotential<0>::set<Lattice, 1>(0.0);
    ChemicalPotential<1>::set<Lattice, 1>(0.0);
    ChemicalPotential<2>::set<Lattice, 1>(0.0);
    Density<>::set<Lattice>(1.0);
    InverseTau<>::set<Lattice>(1.0);

    Pressure<>::set<Lattice>(1);

    Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
        BoundaryLabels<Lattice::NDIM>::template getInstance<Lattice>());
    Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
        OrderParameter<0>::getInstance<Lattice>());
    Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
        OrderParameter<1>::getInstance<Lattice>());
    Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
        Pressure<>::getInstance<Lattice>());

    Lattice::ResetParallelTracking();
    Algorithm lbm(ternary1, ternary2, pressure);
    Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
        LaplacianChemicalPotential<0>::getInstance<Lattice>());
    Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
        LaplacianChemicalPotential<1>::getInstance<Lattice>());
    Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
        LaplacianChemicalPotential<2>::getInstance<Lattice>());
    Lattice::ResetParallelTracking();

    saver.saveHeader(timesteps, saveInterval);
    // for (int i=0;i<lx;i++) std::cout<<
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        TIME = timestep;
        // Save the desired parameters, producing a binary file for each.
        if (timestep % saveInterval == 0) {
            if (mpi.rank == 0) std::cout << "Saving at timestep " << timestep << "." << std::endl;

            saver.saveBoundaries(timestep);
            saver.saveParameter<ChemicalPotential<0>>(timestep, true);
            saver.saveParameter<ChemicalPotential<1>>(timestep, true);
            saver.saveParameter<ChemicalPotential<2>>(timestep, true);
            saver.saveParameter<LaplacianChemicalPotential<0>>(timestep, true);
            saver.saveParameter<LaplacianChemicalPotential<1>>(timestep, true);
            saver.saveParameter<LaplacianChemicalPotential<2>>(timestep, true);
            saver.saveParameter<Pressure<>>(timestep);
            saver.saveParameter<OrderParameter<0>>(timestep, true);
            saver.saveParameter<OrderParameter<1>>(timestep, true);
            saver.saveParameter<Velocity<>, Lattice::NDIM>(timestep);
        }
#ifdef MPIPARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        lbm.evolve();
        if (ternary1.isNan()) {
            std::cout << "here" << std::endl;
            exit(1);
        }
        Lattice::ResetParallelTracking();
        Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
            LaplacianChemicalPotential<0>::getInstance<Lattice>());
        Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
            LaplacianChemicalPotential<1>::getInstance<Lattice>());
        Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
            LaplacianChemicalPotential<2>::getInstance<Lattice>());
        Data_Base<Lattice, typename DefaultTraitPressureTernaryLee<Lattice>::Stencil>::getInstance().communicate(
            Pressure<>::getInstance<Lattice>());
        Lattice::ResetParallelTracking();
    }
}
