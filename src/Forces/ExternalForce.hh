#ifndef FORCES_HEADER
#define FORCES_HEADER
#include "../Parameters.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

class BodyForce{
    public:

        double compute(int xyz,int k) const; //Return force at lattice point k in direction xyz

        void precompute(int k); //Perform any neccessary computations before force is computed

        double computeDensitySource(int k) const; //Calculate any possible source/correction term for density

        double computeVelocitySource(int xyz,int k) const; //Calculate any possible source/correction term for
                                                           //velocity

        void postprocess(int k); //Perform any necessary postprocessing

    private:

        double magnitude=0.001;

        Density<double> m_Density; //Density

};

double BodyForce::compute(int xyz,int k) const{

    return (xyz==0||xyz==2)*magnitude*m_Density.getParameter(k); //Force is just density multiplied by magnitude
                                                                 //in given direction

}

void BodyForce::precompute(int k){ //Not necessary
    
}

void BodyForce::postprocess(int k){ //Not necessary
    
}

double BodyForce::computeDensitySource(int k) const{ //Not necessary

    return 0.0;

}

double BodyForce::computeVelocitySource(int xyz,int k) const{ //Need to correct velocity

    return +compute(xyz,k)*DT/(m_Density.getParameter(0));
    
}

#endif