#include "../../src/Algorithm.hh"
#include "../../src/LBModels/Models.hh"
#include "../../src/Forces/Forces.hh"
#include "../../src/BoundaryModels/Boundaries.hh"
#include "../../src/Data.hh"
#include "../../src/Stencil.hh"
#include "../../src/Global.hh"
#include "../../src/Service.hh"
#include "../../src/Saving.hh"

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
//BASIC MPI
//BINARY
//BOUNCEBACK
//COMMENTS
//CLEANUP (SAVING, EXCEPTIONS etc)
//CLEANUP OLD CODE

//main.cc: This file is just used to run the LBM code and choose how to setup the simulation.

//Modularisation is implemented using trait classes, which contain stencil information, 
//the data type, a tuple of the boundary types and a tuple of forces to be applied in the model.

//Trait class for FlowField Distribution (Navier-Stokes and continuity solver)
struct traitFlowField{
    using Stencil=D2Q9; //Here, D refers to the number of cartesian dimensions
                        //and Q refers to the number of discrete velocity directions.
                        //This naming convention is standard in LBM.
    using Data=Data1<Stencil>; //This will change the "Data" implementation, which will essentially
                               //govern the access of non-local data
    using Boundaries=std::tuple<BounceBack&>; //This will tell the model which boundaries to apply
    using Forces=std::tuple<BodyForce&>; //This will tell the model which forces to apply
};

//Trait class for PhaseField Distribution (Calculates the interface between components)
struct traitPhaseField{
    using Stencil=D2Q9;  
    using Data=Data1<Stencil>;
    using Boundaries=std::tuple<BounceBack&>;
    using Forces=std::tuple<>;
};

int main(){

    data_dir="data/"; //TEMPORARY used to save output
    

    BodyForce o_ForceFlowField; //Create an object of the bodyforce class,
                     //I want people to be able to pass information here through the constructor
    auto ForcesFlowField=GenerateTuple(o_ForceFlowField); //Generate a tuple of all forces passed to this function

    BounceBack o_BoundaryFlowField; //Create an object of the bounceback class
    auto BoundaryFlowField=GenerateTuple(o_BoundaryFlowField); //Generate tuple

    BounceBack o_BoundaryPhaseField; //Create an object of the bounceback class
    auto BoundaryPhaseField=GenerateTuple(o_BoundaryPhaseField); //Generate tuple


    FlowField<traitFlowField> o_FlowField(ForcesFlowField,BoundaryFlowField); //Create an object of the FlowField model class
                                                                      //and pass forces and boundaries to it
    Binary<traitPhaseField> o_PhaseField(BoundaryPhaseField);


    Algorithm<FlowField<traitFlowField>,Binary<traitPhaseField>> LBM(o_FlowField,o_PhaseField);
    

    LBM.initialise(); //Perform necessary initialisation

    for (int timestep=0;timestep<=TIMESTEPS;timestep++){

        LBM.evolve(); //Evolve one timestep

    }
    
    return 0;
}