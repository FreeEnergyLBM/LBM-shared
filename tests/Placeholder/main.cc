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

const int LX=200;
const int LY=200;
const int TIMESTEPS=100;
using Lattice=LatticeProperties<Data1,X_Parallel,LX,LY>;

struct traittt:DefaultTrait<Lattice,FlowField>{

};

int main(int argc, char **argv){

    Lattice l1;

    #ifdef MPIPARALLEL
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NUMPROCESSORS);                              // Store number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &CURPROCESSOR);                              // Store processor IDs
    Parallel<1> initialise(l1);
    #endif
    
    auto Dist1=Model<FlowFieldBinary>(l1);
    auto Dist2=Model<Binary>(l1);

    Algorithm LBM(l1,Dist1,Dist2);
    
    ParameterSave<Density,OrderParameter,Velocity> Saver(l1,"data/");

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
