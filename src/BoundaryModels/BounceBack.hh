#ifndef BOUCNEBACK_HEADER
#define BOUCNEBACK_HEADER
#include "../Parameters.hh"
#include<iostream>

class BounceBack{
    public:
        template <class disttype>
        void compute(disttype& m_Distribution,int k,int idx) const;

        void precompute(int k);

        double computeDensitySource(int k) const;

        double computeVelocitySource(int xyz,int k) const;

        void postprocess(int k);

    private:


};

template <class disttype>
void BounceBack::compute(disttype& m_Distribution,int k,int idx) const{
    m_Distribution.getDistributionPointer(m_Distribution.streamIndex(k,idx))[idx]=m_Distribution.getDistributionPointer(k)[m_Distribution.getOpposite(idx)];
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