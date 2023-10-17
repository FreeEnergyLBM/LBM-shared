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
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (Geometry<Lattice>::getBoundaryType(k) != 4) return;

    const int& normalq = Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        
        
        double cidotnormal = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            cidotnormal += TTraits::Stencil::Ci_xyz(xyz)[idx]*BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection[xyz];
        }
        //if (cidotnormal > 0) TParameter::template get<Lattice>(neighbors[k * Stencil::Q+Stencil::Opposites[idx1]]) = (4*TParameter::template get<Lattice>(k) - TParameter::template get<Lattice>(neighbors[k * Stencil::Q+idx1])) / 3.0;
        
        if(Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == 0 && cidotnormal!=0) {
            
            //distribution.getDistributionPointer(k)[idx] = (4.0*distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] - distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[idx])/3.0;
            distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = (4.0*distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[idx] - distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq), normalq))[idx])/3.0;
            //std::cout<<normalq<<std::endl;
            //distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -distribution.getDistributionPointer(k)[distribution.getOpposite(idx)] + 2*Stencil::Weights[idx]*mInterfaceVal;
        
        }
        else if (cidotnormal==0){
            distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = distribution.getEquilibriumPointer(distribution.streamIndex(k, normalq))[idx];//TTraits::Stencil::Weights[idx]*OrderParameter<>::get<typename TTraits::Lattice>(k);
        }
        //else if (Geometry<Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == 4 ){
        //    distribution.getDistributionPointer(k)[idx] = Stencil::Weights[idx]*(OrderParameter<>::get<Lattice>(k));
        //}

    }    

}

template<class TTraits, class TDistributionType>
inline void ExtrapolationOutflow::communicatePrecompute(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}
