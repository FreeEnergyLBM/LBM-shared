#ifndef BOUCNEBACK_HEADER
#define BOUCNEBACK_HEADER
#include "../Parameters.hh"
#include "../Boundary.hh"
#include<iostream>

class BounceBack{
    public:
        template <class disttype>
        double compute(disttype& m_Distribution,int k) const;

        void precompute(int k);

        double computeDensitySource(int k) const;

        double computeVelocitySource(int xyz,int k) const;

        void postprocess(int k);

    private:


};

template <class disttype>
double BounceBack::compute(disttype& m_Distribution,int k) const{

}

void BounceBack::precompute(int k){
    
}

void BounceBack::postprocess(int k){
    
}

double BounceBack::computeDensitySource(int k) const{

}

double BounceBack::computeVelocitySource(int xyz,int k) const{

}
#endif