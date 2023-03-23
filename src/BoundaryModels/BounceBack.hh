#ifndef BOUCNEBACK_HEADER
#define BOUCNEBACK_HEADER
#include "../Parameters.hh"
#include "../Boundary.hh"
#include<iostream>

class BounceBack{
    public:

        double compute(int xyz,int k) const;

        void precompute(int k);

        double computeDensitySource(int k) const;

        double computeVelocitySource(int xyz,int k) const;

        void postprocess(int k);

    private:

        Distribution<double>& m_Distribution;

};

double BounceBack::compute(int xyz,int k) const{

}

void BodyForce::precompute(int k){
    
}

void BodyForce::postprocess(int k){
    
}

double BodyForce::computeDensitySource(int k) const{

}

double BodyForce::computeVelocitySource(int xyz,int k) const{

}
#endif