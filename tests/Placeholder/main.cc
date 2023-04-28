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

const int LXtemp=100;
const int LYtemp=100;
const int LZtemp=1;
const int TIMESTEPS=10;
enum Dimension {dimension_2D=2,dimension_3D=3};

int main(int argc, char **argv){

    LatticeProperties l1(LXtemp,LYtemp,LZtemp,dimension_2D);

    #ifdef MPIPARALLEL
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NUMPROCESSORS);                              // Store number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &CURPROCESSOR);                              // Store processor IDs
    Parallel<1> initialise(l1);
    #else
    LXdiv=LX
    #endif
    
    Algorithm<FlowFieldBinary<>,Binary<>> LBM(l1);
    
    ParameterSave<Density,OrderParameter,Velocity> Saver(l1,"data/");

    LBM.initialise(); //Perform necessary initialisation
    
    for (int timestep=0;timestep<=TIMESTEPS;timestep++){
        
        if (timestep%50000==0) Saver.Save(timestep);
        LBM.evolve(); //Evolve one timestep

    }
    
    #ifdef MPIPARALLEL
    MPI_Finalize();
    #endif
    
    return 0;

}
