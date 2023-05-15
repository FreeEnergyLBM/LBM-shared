#include <chrono>
#include <iostream>
#ifdef MPIPARALLEL
#include <mpi.h>
#endif
#include "main.hh"
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
//init in x dir
//variable viscosity

/**
 * \file main.cc
 * \brief This file is just used to run the LBM code and choose how to setup the simulation.
 *
 */

//Modularisation is implemented using trait classes, which contain stencil information, 
//the data type, a tuple of the boundary types and a tuple of forces to be applied in the model.

const int LX=100;
const int LY=100;
const int LZ=1;
const int TIMESTEPS=1000;
const int SAVEINTERVAL=100;
using Lattice=LatticeProperties<Data1,X_Parallel,LX,LY,LZ>;

bool fluidLocation(const int k){
    int yy=computeY(LY,LZ,k);
    if ((yy)>LY/2) return true;
    else return false;
}

bool solidLocation(const int k){
    int yAtCurrentk = computeY(LY, LZ, k);
    if (yAtCurrentk <= 1 || yAtCurrentk >= LY - 2 ) return true;
    else return false;
}

int main(int argc, char **argv){
    
    #ifdef MPIPARALLEL
    int provided;
    MPI_Init_thread(&argc, &argv,MPI_THREAD_FUNNELED,&provided);
    MPI_Comm_size(MPI_COMM_WORLD, &NUMPROCESSORS);                              // Store number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &CURPROCESSOR);                              // Store processor IDs
    
    Parallel<Lattice,1> initialise;
    #endif

    FlowFieldBinary<Lattice> Dist1;
    Binary<Lattice> Dist2;

    Dist2.setTau1(0.51);

    Dist1.getForce<BodyForce>().setMagnitudeX(0.000001);

    OrderParameter<Lattice> orderparam;
    orderparam.set(fluidLocation,-1.0);

    SolidLabels<Lattice> solid;
    solid.set(solidLocation,true);

    Algorithm LBM(Dist1,Dist2);

    ParameterSave<Lattice,Density,OrderParameter,Velocity> Saver("data/");
    Saver.SaveHeader(TIMESTEPS,SAVEINTERVAL);
    
    LBM.initialise(); //Perform necessary initialisation
    
    auto t0=std::chrono::system_clock::now();
    
    #pragma omp parallel
    {
    for (int timestep=0;timestep<=TIMESTEPS;timestep++){
        #pragma omp master
        {
        if (timestep%100==0) Saver.Save(timestep);
        }
        LBM.evolve(); //Evolve one timestep
        
    }
    }
    auto tend=std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds=tend-t0; 

    if(CURPROCESSOR==0)std::cout<<"RUNTIME: "<<elapsed_seconds.count()<<" "<<LX<<" "<<LY<<" "<<std::endl;

    #ifdef MPIPARALLEL
    MPI_Finalize();
    #endif
    
    return 0;

}
