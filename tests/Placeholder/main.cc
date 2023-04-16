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

//TODO BEFORE HACKATHON
//EXCEPTIONS & CLEANUP
//CLEANUP OLD CODE

/**
 * \file main.cc
 * \brief This file is just used to run the LBM code and choose how to setup the simulation.
 *
 */

//Modularisation is implemented using trait classes, which contain stencil information, 
//the data type, a tuple of the boundary types and a tuple of forces to be applied in the model.

//Trait class for FlowField Distribution (Navier-Stokes and continuity solver)

struct traitFlowField{
    using Stencil=D2Q9; //Here, D refers to the number of cartesian dimensions
                        //and Q refers to the number of discrete velocity directions.
                        //This naming convention is standard in LBM.
    #ifdef MPIPARALLEL
    using Parallel=X_Parallel<Stencil,NO_NEIGHBOR>;
    #else
    using Parallel=No_Parallel;
    #endif
    using Data=Data1<Stencil,Parallel>; //This will change the "Data" implementation, which will essentially
                               //govern the access of non-local data
    using Boundaries=std::tuple<BounceBack>; //This will tell the model which boundaries to apply
    using Forces=std::tuple<BodyForce,ChemicalForce>; //This will tell the model which forces to apply
};


//Trait class for PhaseField Distribution (Calculates the interface between components)
struct traitPhaseField{
    using Stencil=D2Q9;  
    #ifdef MPIPARALLEL
    using Parallel=X_Parallel<Stencil,NO_NEIGHBOR>;
    #else
    using Parallel=No_Parallel;
    #endif
    using Data=Data1<Stencil,Parallel>;
    using Boundaries=std::tuple<BounceBack>;
    using Forces=std::tuple<OrderParameterGradients<CentralXYZ<Stencil,Parallel>>>;
};

int main(int argc, char **argv){

    #ifdef MPIPARALLEL
    //#ifdef OMPPARALLEL
    //MPI_Init_thread(&argc, &argv,MPI_THREAD_MULTIPLE, &prov);                                            // Initialise parallelisation based on arguments given
    //#else
    MPI_Init(&argc, &argv);
    //#endif
    MPI_Comm_size(MPI_COMM_WORLD, &NUMPROCESSORS);                              // Store number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &CURPROCESSOR);                              // Store processor IDs
    Parallel<NO_NEIGHBOR> initialise;
    #ifdef OMPPARALLEL
    #pragma omp parallel
    #pragma omp critical
    {
        #pragma omp master
        if(CURPROCESSOR==0) std::cout<<"MPI Processes: "<<NUMPROCESSORS<<" Threads: "<<omp_get_num_threads()<<std::endl;
    
    }
    #endif
    #endif
    
    system("mkdir data");
    DATA_DIR="data/"; //TEMPORARY used to save output

    FlowFieldBinary<traitFlowField> o_FlowField; //Create an object of the FlowField model class
                                                                     //and pass forces and boundaries to it
    Binary<traitPhaseField> o_PhaseField;

    Algorithm<FlowFieldBinary<traitFlowField>,Binary<traitPhaseField>> LBM(o_FlowField,o_PhaseField);
    
    LBM.initialise(); //Perform necessary initialisation
    auto t0=std::chrono::system_clock::now();
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
    auto tend=std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds=tend-t0; 
    std::cout<<"RUNTIME: "<<elapsed_seconds.count()<<" "<<LX<<" "<<LY<<std::endl;
    
    #ifdef MPIPARALLEL
    MPI_Finalize();
    #endif
    #ifdef OMPPARALLEL
    //std::cout<<"RUNTIME: "<<TOTALTIME<<std::endl;
    #endif
    
    return 0;
}