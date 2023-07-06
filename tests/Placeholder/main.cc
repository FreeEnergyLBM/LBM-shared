#include<../../src/lbm.hh>
#include <chrono>
#include <iostream>
#include <math.h>

//TODO
//Swapping based streaming
//BLAS
//non-cuboid grid - sparse matrices
//buffer size
//adative mesh
//floats

/**
 * \file main.cc
 * \brief This file is just used to run the LBM code and choose how to setup the simulation.
 *
 */

//Set up the lattice, including the resolution and data/parallelisation method
const int LX = 100; //Size of domain in x direction
const int LY = 100; //Size of domain in y direction
const int LZ = 1; //Size of domain in z direction (Can also not specify LZ if it is 1)
using Lattice = LatticeProperties<Data1, X_Parallel<1>, LX, LY, LZ>;

const int TIMESTEPS =1000; //Number of iterations to perform
const int SAVEINTERVAL = 1000; //Interval to save global data

const double RADIUS=30.; //Droplet radius

double A = 0.00025;
double kappa = 0.0004;

double tau1 = 1.;
double tau2 = 0.6;

//User defined function to define some fluid initialisation (optional)
bool fluidLocation(int k) {

    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(LY, LZ, k);

    double rr2 = (xx - LX/2.) * (xx - LX/2.) + (yy - LX/2.) * (yy - LX/2.);

    if (rr2 < RADIUS*RADIUS) return true;
    else return false;
    
}

double fluidLocationSmooth(int k) {

    int yy = computeY(LY, LZ, k);

    return tanh((yy-LY/2.)/(sqrt(2*kappa/A)));
}

double fluidLocationSmoothDroplet(int k) {

    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(LY, LZ, k);

    double rr2 = (xx - LX/2.) * (xx - LX/2.) + (yy - 0*LX/2.) * (yy - 0*LX/2.);
    return tanh((sqrt(rr2)-RADIUS)/(sqrt(2*kappa/A)));

}

//User defined function to define some solid initialisation
bool solidLocation(const int k) {

    int yAtCurrentk = computeY(LY, LZ, k);

    if (yAtCurrentk <= 1 || yAtCurrentk >= LY - 2) return true;
    else return false;

}

using flowfieldtraits = DefaultTraitFlowFieldBinary<Lattice> :: SetCollisionModel<MRT>;

int main(int argc, char **argv){
    mpi.init();

    //Chosen models
    FlowFieldBinary<Lattice,flowfieldtraits> FlowFieldModel;
    Binary<Lattice> ComponentSeparationModel;

    FlowFieldModel.setTau1(tau1);
    FlowFieldModel.setTau2(tau2);

    ComponentSeparationModel.setTau1(tau1);
    ComponentSeparationModel.setTau2(tau2);

    ComponentSeparationModel.getPreProcessor<ChemicalPotentialCalculatorBinary>().setA(A);
    ComponentSeparationModel.getPreProcessor<ChemicalPotentialCalculatorBinary>().setKappa(kappa);

    //ComponentSeparationModel.getPreProcessor<LinearWetting>().setPrefactor(A,kappa);
    ComponentSeparationModel.getPreProcessor<CubicWetting>().setThetaDegrees(45);

    FlowFieldModel.getForce<BodyForce<>>().setMagnitudeX(0.0000001);

    OrderParameter<>::set<Lattice>(fluidLocationSmooth);
    
    SolidLabels<>::set<Lattice>(solidLocation,1); //Set solid to true where the function we defined previously is true (false by default so don't need to specify this)

    //Algorithm that will combine the models and run them in order
    Algorithm LBM(FlowFieldModel,ComponentSeparationModel); //Create LBM object with the two models we have initialised

    //Saving class
    ParameterSave<Lattice> Saver("data/");
    Saver.SaveHeader(TIMESTEPS,SAVEINTERVAL); //Save header with lattice information (LX, LY, LZ, NDIM (2D or 3D), TIMESTEPS, SAVEINTERVAL)
    
    LBM.initialise(); //Perform necessary initialisation for the models in LBM

    //Loop over timesteps
    for (int timestep=0;timestep<=TIMESTEPS;timestep++) {

        if (timestep%SAVEINTERVAL==0) {
            std::cout<<"Saving at timestep "<<timestep<<"."<<std::endl;
            Saver.SaveParameter<OrderParameter<>>(timestep);
            Saver.SaveParameter<Velocity<>,Lattice::NDIM>(timestep);
        }
         
        LBM.evolve(); //Evolve one timestep of the algorithm
        
    }

    return 0;

}
