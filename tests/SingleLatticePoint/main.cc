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

#include <iostream>
#include <omp.h>
#include <mpi.h>

//TODO BEFORE HACKATHON
//EXCEPTIONS & CLEANUP
//CLEANUP OLD CODE

//main.cc: This file is just used to run the LBM code and choose how to setup the simulation.

//Modularisation is implemented using trait classes, which contain stencil information, 
//the data type, a tuple of the boundary types and a tuple of forces to be applied in the model.

//Trait class for FlowField Distribution (Navier-Stokes and continuity solver)

struct traitFlowField{
    using Stencil=D2Q9; //Here, D refers to the number of cartesian dimensions
                        //and Q refers to the number of discrete velocity directions.
                        //This naming convention is standard in LBM.
    using Parallel=X_Parallel<Stencil,NO_NEIGHBOR>;
    using Data=Data1<Stencil,Parallel>; //This will change the "Data" implementation, which will essentially
                               //govern the access of non-local data
    using Boundaries=std::tuple<BounceBack>; //This will tell the model which boundaries to apply
    //using Forces=std::tuple<BodyForce,ChemicalForce>; //This will tell the model which forces to apply
    using Forces=std::tuple<>;
};


//Trait class for PhaseField Distribution (Calculates the interface between components)
struct traitPhaseField{
    using Stencil=D2Q9;  
    using Parallel=X_Parallel<Stencil,NO_NEIGHBOR>;
    using Data=Data1<Stencil,Parallel>;
    using Boundaries=std::tuple<BounceBack>;
    //using Forces=std::tuple<OrderParameterGradients<CentralXYZ<Stencil,Parallel>>>;
    using Forces=std::tuple<>;
};



int main(int argc, char **argv){
    
    #ifdef MPIPARALLEL
    MPI_Init(&argc, &argv);                                            // Initialise parallelisation based on arguments given
    MPI_Comm_size(MPI_COMM_WORLD, &NUMPROCESSORS);                              // Store number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &CURPROCESSOR);                              // Store processor IDs
    Parallel<NO_NEIGHBOR> initialise;
    #endif

    system("mkdir data");
    DATA_DIR="data/"; //TEMPORARY used to save output

    FlowFieldBinary<traitFlowField> o_FlowField; //Create an object of the FlowField model class
                                                                     //and pass forces and boundaries to it
    Binary<traitPhaseField> o_PhaseField;

    Algorithm<FlowFieldBinary<traitFlowField>,Binary<traitPhaseField>> LBM(o_FlowField,o_PhaseField);
    LBM.initialise(); //Perform necessary initialisation
    
    for (int timestep=0;timestep<=TIMESTEPS;timestep++){

        LBM.evolve(); //Evolve one timestep

        if (timestep%SAVEINTERVAL==0) {
            if(CURPROCESSOR==0) std::cout<<"SAVING at timestep "<<timestep<<""<<std::endl;
            Density<double>::save("density",timestep);
            OrderParameter<double>::save("orderparameter",timestep);
            ChemicalPotential<double>::save("chemicalpotential",timestep);
            Velocity<double,NDIM>::save("velocity",timestep);
        }
        
    }
    
    #ifdef MPIPARALLEL
    MPI_Finalize();
    #endif
    #ifdef OMPPARALLEL
    std::cout<<TOTALTIME<<std::endl;
    #endif
    return 0;
}