#include "../../src/Algorithm.hh"
#include "../../src/LBModels/Models.hh"
#include "../../src/Forces/Forces.hh"
#include "../../src/BoundaryModels/Boundaries.hh"
#include "../../src/GradientStencils/GradientStencils.hh"
#include "../../src/Data.hh"
#include "../../src/Stencil.hh"
#include "../../src/Global.hh"
#include "../../src/Service.hh"
#include "../../src/Saving.hh"
#include "../../src/Parallel.hh"
#include "../../src/Parameters.hh"

#include <iostream>

//Code wishlist:

//	Grid refinement
//	OpenMP
//	Advanced data - think about affinity
//	DDF shifting
//  Use BLAS
// Ci script
// Debugger
// Sphinx
// resize arrays to not include solid
// more things at compile time

//TODO BEFORE HACKATHON
//BINARY
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
    using Data=Data1<Stencil,X_Parallel<Stencil,1>>; //This will change the "Data" implementation, which will essentially
                               //govern the access of non-local data
    using Boundaries=std::tuple<BounceBack&>; //This will tell the model which boundaries to apply
    using Forces=std::tuple<BodyForce&,ChemicalForce&>; //This will tell the model which forces to apply
};

//Trait class for PhaseField Distribution (Calculates the interface between components)
struct traitPhaseField{
    using Stencil=D2Q9;  
    using Data=Data1<Stencil,X_Parallel<Stencil,1>>;
    using Boundaries=std::tuple<BounceBack&>;
    using Forces=std::tuple<OrderParameterGradients<CentralXYZ<Data_Base<Stencil,X_Parallel<Stencil,1>>>>&>;
};

int main(int argc, char **argv){
    
    #ifdef PARALLEL
    MPI_Init(&argc, &argv);                                            // Initialise parallelisation based on arguments given
    MPI_Comm_size(MPI_COMM_WORLD, &NUMPROCESSORS);                              // Store number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &CURPROCESSOR);                              // Store processor IDs
    Parallel<1> initialise;
    #endif
    
    DATA_DIR="data/"; //TEMPORARY used to save output
    BounceBack o_BoundaryFlowField; //Create an object of the bounceback class
    auto BoundaryFlowField=std::tie(o_BoundaryFlowField); //Generate tuple

    BounceBack o_BoundaryPhaseField; //Create an object of the bounceback class
    auto BoundaryPhaseField=std::tie(o_BoundaryPhaseField); //Generate tuple
    BodyForce o_ForceFlowField; //Create an object of the bodyforce class,
                     //I want people to be able to pass information here through the constructor
    ChemicalForce o_ForceFlowFieldChemical;
    auto ForcesFlowField=std::tie(o_ForceFlowField,o_ForceFlowFieldChemical); //Generate a tuple of all forces passed to this function

    OrderParameterGradients<CentralXYZ<Data_Base<traitPhaseField::Stencil,X_Parallel<traitPhaseField::Stencil,1>>>> o_BinaryGradients;
    auto ForcesPhaseField=std::tie(o_BinaryGradients);




    FlowFieldBinary<traitFlowField> o_FlowField(ForcesFlowField,BoundaryFlowField); //Create an object of the FlowField model class
                                                                     //and pass forces and boundaries to it
    Binary<traitPhaseField> o_PhaseField(ForcesPhaseField,BoundaryPhaseField);

     
    Algorithm<FlowFieldBinary<traitFlowField>,Binary<traitPhaseField>> LBM(o_FlowField,o_PhaseField);

    LBM.initialise(); //Perform necessary initialisation
    
    for (int timestep=0;timestep<=TIMESTEPS;timestep++){
        
        LBM.evolve(); //Evolve one timestep
        
        if (timestep%10000==0) {
            Density<double>::save("density",timestep);
            OrderParameter<double>::save("orderparameter",timestep);
            ChemicalPotential<double>::save("chemicalpotential",timestep);
            Velocity<double,NDIM>::save("velocity",timestep);
        }
    }
    
    #ifdef PARALLEL
    MPI_Finalize();
    #endif
    
    return 0;
}