#include "../../src/lbm.hh"
#include <chrono>
#include <iostream>
#ifdef OMPPARALLEL
#include <omp.h>
#endif
#ifdef MPIPARALLEL
#include <mpi.h>
#endif

//TODO
//Swapping based streaming
//BLAS
//non-cuboid grid - sparse matrices
//buffer size
//adative mesh
//floats
//MPI object - get rid of ifdefs
//Sort out globals

/**
 * \file main.cc
 * \brief This file is just used to run the LBM code and choose how to setup the simulation.
 *
 */

//Modularisation is implemented using trait classes, which contain stencil information, 
//the data type, a tuple of the boundary types and a tuple of forces to be applied in the model.

const int LX=100;
const int LY=100;
const int TIMESTEPS=10000;
using Lattice=LatticeProperties<Data1,X_Parallel,LX,LY>;
//using Lattice=LatticePropertiesRuntime<Data1,X_Parallel,2>;

auto& GETPROPERTIES(){
    return getGlobal<Lattice>();
}

int main(int argc, char **argv){

    #ifdef MPIPARALLEL
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NUMPROCESSORS);                              // Store number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &CURPROCESSOR);                              // Store processor IDs
    Parallel<1> initialise;
    #endif

    FlowFieldBinary Dist1;
    Binary Dist2;

    Dist1.getForce<BodyForce>().setMagnitude(0.0001);

    Algorithm LBM(Dist1,Dist2);
    
    ParameterSave<Density,OrderParameter,Velocity> Saver("data/");

    LBM.initialise(); //Perform necessary initialisation
    auto t0=std::chrono::system_clock::now();

    for (int timestep=0;timestep<=TIMESTEPS;timestep++){
        
        if (timestep%50000==0) Saver.Save(timestep);
        LBM.evolve(); //Evolve one timestep

    }

    auto tend=std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds=tend-t0; 
    std::cout<<"RUNTIME: "<<elapsed_seconds.count()<<" "<<LX<<" "<<LY<<std::endl;

    #ifdef MPIPARALLEL
    MPI_Finalize();
    #endif
    
    return 0;

}
