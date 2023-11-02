#pragma once
#include "BoundaryBase.hh"
#include "../Geometry.hh"
#include<iostream>


class BounceBack : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

    private:

};

template<class TTraits, class TDistributionType>
inline void BounceBack::compute(TDistributionType& distribution, int k) {

    if (Geometry<typename TTraits::Lattice>::getBoundaryType(k) != 1) return;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

        if(Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == 0 ) {
            
            distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx));
        
        }

    }    

}