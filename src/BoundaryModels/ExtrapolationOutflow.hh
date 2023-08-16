#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class ExtrapolationOutflow : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicatePrecompute(){};

        template<class TTraits, class TDistributionType>
        inline void communicatePrecompute(TDistributionType& mDistribution);

    private:

};

template<class TTraits, class TDistributionType>
inline void ExtrapolationOutflow::compute(TDistributionType& distribution, int k) { //CHANGE THIS SO YOU DONT NEED TO COMMUNICATE

    if (Geometry<typename TTraits::Lattice>::getBoundaryType(k) != 4) return;

    const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<>::get<typename TTraits::Lattice>(k).NormalDirection)->second;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

        if(Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == 0 ) {
            
            distribution.getDistributionPointer(k)[idx] = 4.0/3.0*distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] - distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[idx]/3.0;
            //distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -distribution.getDistributionPointer(k)[distribution.getOpposite(idx)] + 2*TTraits::Stencil::Weights[idx]*mInterfaceVal;
        
        }

    }    

}

template<class TTraits, class TDistributionType>
inline void ExtrapolationOutflow::communicatePrecompute(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}