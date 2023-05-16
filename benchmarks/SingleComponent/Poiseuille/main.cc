#include<../../../src/lbm.hh>
#include <chrono>
#include <iostream>
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
//Adjust parallelisation based on solid
//template function to add forces

/**
 * \file main.cc
 * \brief This file is just used to run the LBM code and choose how to setup the simulation.
 *
 */

//Set up the lattice, including the resolution and data/parallelisation method
const int LX = 10; //Size of domain in x direction
const int LY = 100; //Size of domain in y direction
using Lattice = LatticeProperties<Data1, X_Parallel, LX, LY>;

const int TIMESTEPS = 10000; //Number of iterations to perform
const int SAVEINTERVAL = 10000; //Interval to save global data

//User defined function to define some fluid initialisation (optional)

//User defined function to define some solid initialisation
bool solidLocation(const int k) {

    int yAtCurrentk = computeY(LY, 1, k);

    if (yAtCurrentk <= 1 || yAtCurrentk >= LY - 2) return true;
    else return false;

}


//Traits of the model. We are using the defaults but adding a body force.
using PoiseuilleTrait = DefaultTraitBinary<Lattice> ::AddAddOn<BodyForce<Lattice>>;

int main(int argc, char **argv){
    
    #ifdef MPIPARALLEL
    //MPI initialisation
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &NUMPROCESSORS);                              // Store number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &CURPROCESSOR);                              // Store processor IDs
    
    Parallel<Lattice,1> initialise; //Initialise parallelisation
    #endif

    //Chosen models
    FlowField<Lattice,PoiseuilleTrait> Model; //Flowfield

    Model.getAddOn<BodyForce<Lattice>>().setMagnitudeX(0.000001); //Get object of body force and then set the magnitude

    SolidLabels<Lattice> solid;
    solid.set(solidLocation,true); //Set solid to true where the function we defined previously is true (false by default so don't need to specify this)

    //Algorithm that will combine the models and run them in order
    Algorithm LBM(Model); //Create LBM object with the model we have initialised

    //Saving class
    ParameterSave<Lattice,Density,OrderParameter,Velocity> Saver("data/"); //Specify the lattice, the parameters you want to save and the  data directory
    Saver.SaveHeader(TIMESTEPS,SAVEINTERVAL); //Save header with lattice information (LX, LY, LZ, NDIM (2D or 3D), TIMESTEPS, SAVEINTERVAL)
    
    LBM.initialise(); //Perform necessary initialisation for the models in LBM
    
    //Loop over timesteps
    for (int timestep=0;timestep<=TIMESTEPS;timestep++) {

        if (timestep%SAVEINTERVAL==0) Saver.Save(timestep);

        LBM.evolve(); //Evolve one timestep of the algorithm
        
    }

    #ifdef MPIPARALLEL
    MPI_Finalize(); //MPI finalisation
    #endif
    
    return 0;

}
