#ifndef FORCES_HEADER
#define FORCES_HEADER
#include "../Parameters.hh"
#include<iostream>

class BodyForce{
    public:

        double compute(int xyz,int k) const;

        void precompute(int k);

        double computeDensitySource(int k) const;

        double computeVelocitySource(int xyz,int k) const;

        void postprocess(int k);

    private:

        double magnitude=0.001;

        Density<double> m_Density;

};

double BodyForce::compute(int xyz,int k) const{

    return (xyz==0||xyz==2)*magnitude*m_Density.getParameter(k);

}

void BodyForce::precompute(int k){
    
}

void BodyForce::postprocess(int k){
    
}

double BodyForce::computeDensitySource(int k) const{

    return 0.0;

}

double BodyForce::computeVelocitySource(int xyz,int k) const{

    return +compute(xyz,k)*DT/(m_Density.getParameter(0));
    
}

#endif