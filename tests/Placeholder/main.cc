#include<../../src/lbm.hh>
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
//Adjust parallelisation basedo n solid
//template function to add forces

/**
 * \file main.cc
 * \brief This file is just used to run the LBM code and choose how to setup the simulation.
 *
 */

//Set up the lattice, including the resolution and data/parallelisation method
const int LX = 100; //Size of domain in x direction
const int LY = 100; //Size of domain in y direction
const int LZ = 1; //Size of domain in z direction (Can also not specify LZ if it is 1)
using Lattice = LatticeProperties<Data1, X_Parallel, LX, LY, LZ>;

const int TIMESTEPS = 1000; //Number of iterations to perform
const int SAVEINTERVAL = 100; //Interval to save global data

//User defined function to define some fluid initialisation (optional)
bool fluidLocation(const int k) {

    int yy = computeY(LY, LZ, k);
    
    if (yy > LY / 2) return true;
    else return false;
    
}

//User defined function to define some solid initialisation
bool solidLocation(const int k) {

    int yAtCurrentk = computeY(LY, LZ, k);

    if (yAtCurrentk <= 1 || yAtCurrentk >= LY - 2) return true;
    else return false;

}

int main(int argc, char **argv){
    
    #ifdef MPIPARALLEL
    //MPI initialisation
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &NUMPROCESSORS);                              // Store number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &CURPROCESSOR);                              // Store processor IDs
    
    Parallel<Lattice,1> initialise; //Initialise parallelisation
    #endif

    //Chosen models
    FlowFieldBinary<Lattice> Model1; //Flowfield (navier stokes solver) that can be used with the binary model (there are nuances with this model)
    Binary<Lattice> Model2; //Binary model with hybrid equilibrium and forcing term

    Model1.getForce<BodyForce>().setMagnitudeX(0.000001); //Get object of body force and then set the magnitude
    Model2.setTau1(0.51); //Set the relaxation time of fluid 1

    OrderParameter<Lattice> orderparam;
    orderparam.set(fluidLocation, -1.0, 1.0); //Set fluid to -1 where the function we defined previously is true and 1.0 where it is false

    SolidLabels<Lattice> solid;
    solid.set(solidLocation,true); //Set solid to true where the function we defined previously is true (false by default so don't need to specify this)

    //Algorithm that will combine the models and run them in order
    Algorithm LBM(Model1,Model2); //Create LBM object with the two models we have initialised

    //Saving class
    ParameterSave<Lattice,Density,OrderParameter,Velocity> Saver("data/"); //Specify the lattice, the parameters you want to save and the  data directory
    Saver.SaveHeader(TIMESTEPS,SAVEINTERVAL); //Save header with lattice information (LX, LY, LZ, NDIM (2D or 3D), TIMESTEPS, SAVEINTERVAL)
    
    LBM.initialise(); //Perform necessary initialisation for the models in LBM
    
    auto t0=std::chrono::system_clock::now(); //Start a timer
    
    //Loop over timesteps
    for (int timestep=0;timestep<=TIMESTEPS;timestep++) {

        if (timestep%100==0) Saver.Save(timestep);

        LBM.evolve(); //Evolve one timestep of the algorithm
        
    }

    auto tend=std::chrono::system_clock::now(); //End timer
    std::chrono::duration<double> elapsed_seconds=tend-t0; //Work out total time (end minus start)
    if(CURPROCESSOR==0)std::cout<<"RUNTIME: "<<elapsed_seconds.count()<<" "<<LX<<" "<<LY<<" "<<std::endl; //Print runtime

    #ifdef MPIPARALLEL
    MPI_Finalize(); //MPI finalisation
    #endif
    
    return 0;

}
