#include "../../src/Algorithm.hh"
#include "../../src/LBModels/Models.hh"
#include "../../src/Forces/Forces.hh"
#include "../../src/BoundaryModels/Boundaries.hh"
#include "../../src/GradientStencils/GradientStencils.hh"
#include "../../src/Data.hh"
#include "../../src/Stencil.hh"
#include "../../src/Global.hh"
#include "../../src/Service.hh"
#include "../../src/Parallel.hh"
#include "../../src/Parameters.hh"
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

int main(int argc, char **argv){

    #ifdef MPIPARALLEL
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NUMPROCESSORS);                              // Store number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &CURPROCESSOR);                              // Store processor IDs
    Parallel<NO_NEIGHBOR> initialise;
    #endif
    
    system("mkdir data");
    DATA_DIR="data/"; //TEMPORARY used to save output

    Algorithm<FlowFieldBinary<>,Binary<>> LBM;

    LBM.initialise(); //Perform necessary initialisation

    for (int timestep=0;timestep<=TIMESTEPS;timestep++){

        if (timestep%SAVEINTERVAL==0) {
            if(CURPROCESSOR==0) std::cout<<"SAVING at timestep "<<timestep<<""<<std::endl;
            Density<double>::save("density",timestep);
            OrderParameter<double>::save("orderparameter",timestep);
            ChemicalPotential<double>::save("chemicalpotential",timestep);
            Velocity<double,NDIM>::save("velocity",timestep);
        }

        LBM.evolve(); //Evolve one timestep

    }
    
    #ifdef MPIPARALLEL
    MPI_Finalize();
    #endif
    #ifdef OMPPARALLEL
    //std::cout<<"RUNTIME: "<<TOTALTIME<<std::endl;
    #endif
    
    return 0;
}