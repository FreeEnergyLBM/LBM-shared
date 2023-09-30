#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class Dirichlet : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicatePrecompute(){};

        template<class TTraits, class TDistributionType>
        inline void communicatePrecompute(TDistributionType& mDistribution);

        inline void setInterfaceVal(double val) {mInterfaceVal=val;};

        inline void setInterfaceID(int id) {mInterfaceID=id;};

    private:

        double mInterfaceVal;
        int mInterfaceID;

};

template<class TTraits, class TDistributionType>
inline void Dirichlet::compute(TDistributionType& distribution, int k) { //CHANGE THIS SO YOU DONT NEED TO COMMUNICATE

    if (Geometry<typename TTraits::Lattice>::getBoundaryType(k) != mInterfaceID) return;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

        if(Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == 0 ) {
            
            distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -distribution.getDistributionOldPointer(distribution.streamIndex(k, idx))[distribution.getOpposite(idx)] + 2*TTraits::Stencil::Weights[idx]*mInterfaceVal;
            //distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -distribution.getDistributionPointer(k)[distribution.getOpposite(idx)] + 2*TTraits::Stencil::Weights[idx]*mInterfaceVal;
        
        }

    }    

}

template<class TTraits, class TDistributionType>
inline void Dirichlet::communicatePrecompute(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}