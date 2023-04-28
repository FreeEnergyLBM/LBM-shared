#ifndef EXFORCE_HEADER
#define EXFORCE_HEADER
#include "../Parameters.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

class BodyForce{
    public:

        BodyForce(LatticeProperties& properties):m_Density(properties),m_Properties(properties){}

        double compute(int xyz,int k) const; //Return force at lattice point k in direction xyz

        void precompute(int k); //Perform any neccessary computations before force is computed

        double computeDensitySource(int k) const; //Calculate any possible source/correction term for density

        double computeVelocitySource(int xyz,int k) const; //Calculate any possible source/correction term for
                                                           //velocity

        void postprocess(int k); //Perform any necessary postprocessing

    private:

        LatticeProperties& m_Properties;
        const double& DT=m_Properties.m_DT;

        double magnitude=0.00000001;//1;

        Density m_Density; //Density

};

double BodyForce::compute(int xyz,int k) const{

    return (xyz==0)*magnitude*m_Density.getParameter(k); //Force is just density multiplied by magnitude
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

    return +compute(xyz,k)*DT/(2.0*m_Density.getParameter(k));
    
}

#endif