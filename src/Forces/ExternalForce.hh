#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

class BodyForce{
    public:
        template<template<class,class> class data,template<class,int> class parallel,int lx, int ly,int lz=1>
        BodyForce(LatticeProperties<data,parallel,lx,ly,lz>& properties)
            : DT(properties.m_DT),
              m_Density(properties){}

        BodyForce(const BodyForce& other)
            : DT(other.DT),
              m_Density(other.m_Density){}

        double compute(int xyz,int k) const; //Return force at lattice point k in direction xyz

        void precompute(int k); //Perform any neccessary computations before force is computed

        double computeDensitySource(int k) const; //Calculate any possible source/correction term for density

        double computeVelocitySource(int xyz,int k) const; //Calculate any possible source/correction term for
                                                           //velocity

        void postprocess(int k); //Perform any necessary postprocessing

        const double& DT;

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