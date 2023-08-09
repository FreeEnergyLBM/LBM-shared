#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class Dirichlet : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        inline void setInterfaceVal(double val) {mInterfaceVal=val;};

        inline void setInterfaceID(int id) {mInterfaceID=id;};

    private:

        double mInterfaceVal;
        int mInterfaceID;

};

template<class TTraits, class TDistributionType>
inline void Dirichlet::compute(TDistributionType& distribution, int k) {

    if (Geometry<typename TTraits::Lattice>::getBoundaryType(k) != mInterfaceID) return;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

        if(Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == 0 ) {
            
            distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -distribution.getDistributionPointer(k)[distribution.getOpposite(idx)] + 2*TTraits::Stencil::Weights[idx]*mInterfaceVal;
        
        }

    }    

}