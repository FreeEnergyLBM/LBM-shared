#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class Temp0Outflow : public BoundaryBase {
    public:

        Temp0Outflow() : BoundaryBase(4) {}

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);


};

template<class TTraits, class TDistributionType>
inline void Temp0Outflow::compute(TDistributionType& distribution, int k) { //CHANGE THIS SO YOU DONT NEED TO COMMUNICATE
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    const int& normalq = 2;//Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;

    for (int idx = 0; idx < Stencil::Q; idx++) {
        
        //if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;
        if (BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(distribution.streamIndex(k, idx)).Id!=0) continue;

        distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = Stencil::Weights[idx] * Density<>::get<Lattice>(distribution.streamIndex(k, normalq));
        //#pragma omp critical
        //{
        //std::cout<<k<<" "<<(distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[idx])<<std::endl;
        //}
        

    }    

}

template<class TTraits, class TDistributionType>
inline void Temp0Outflow::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
    Lattice::communicateDistributionAllOld(distribution);

}
