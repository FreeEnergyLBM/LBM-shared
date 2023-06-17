#include<../../../src/lbm.hh>
#include <chrono>
#include <iostream>

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
const int LX = 100; //Size of domain in x direction
const int LY = 100; //Size of domain in y direction
using Lattice = LatticeProperties<Data1, X_Parallel<1>, LX, LY>;

const int TIMESTEPS = 10000; //Number of iterations to perform
const int SAVEINTERVAL = 10000; //Interval to save global data

const int RADIUS=20; //Droplet radius

//User defined function to define some fluid initialisation (optional)
bool fluidLocation(const int k) {

    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(LY, 1, k);

    const int rr2 = (xx - LX / 2) * (xx - LX / 2) + (yy - LY / 2) * (yy - LY / 2);

    if (rr2 < RADIUS*RADIUS) return true;
    else return false;
    
}

int main(int argc, char **argv){
    mpi.init();

    //Chosen models
    FlowFieldBinary<Lattice> Model1; //Flowfield (navier stokes solver) that can be used with the binary model (there are nuances with this model)
    Binary<Lattice> Model2; //Binary model with hybrid equilibrium and forcing term

    //Set parameters relating to the interface width and surface tension
    Model1.getAddOn<ChemicalForce<Lattice>>().setA(0.00015);
    Model1.getAddOn<ChemicalForce<Lattice>>().setKappa(0.0003);

    OrderParameter<Lattice> orderparam;
    orderparam.set(fluidLocation, -1.0, 1.0); //Set fluid to -1 where the function we defined previously is true and 1.0 where it is false

    //Algorithm that will combine the models and run them in order
    Algorithm LBM(Model1,Model2); //Create LBM object with the two models we have initialised

    //Saving class
    ParameterSave<Lattice,Density,OrderParameter,Velocity> Saver("data/"); //Specify the lattice, the parameters you want to save and the  data directory
    Saver.SaveHeader(TIMESTEPS,SAVEINTERVAL); //Save header with lattice information (LX, LY, LZ, NDIM (2D or 3D), TIMESTEPS, SAVEINTERVAL)
    
    LBM.initialise(); //Perform necessary initialisation for the models in LBM
    
    //Loop over timesteps
    for (int timestep=0;timestep<=TIMESTEPS;timestep++) {

        if (timestep%SAVEINTERVAL==0) Saver.Save(timestep);

        LBM.evolve(); //Evolve one timestep of the algorithm
        
    }
    
    return 0;

}
