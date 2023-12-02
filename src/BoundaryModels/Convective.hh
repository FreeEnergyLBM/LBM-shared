#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class Convective : public BoundaryBase {
    public:

        Convective() { this->setInterfaceID(4); }

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);


};

template<class TTraits, class TDistributionType>
inline void Convective::compute(TDistributionType& distribution, int k) { //CHANGE THIS SO YOU DONT NEED TO COMMUNICATE
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    const int& normalq = Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;
    const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;
    //#pragma omp critical
    //{
    //std::cout<<normalq<<std::endl;
    //}
    //std::cout<<normalq<<std::endl;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        
        //if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;
        if (BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(distribution.streamIndex(k, idx)).Id!=0) continue;

        double normalvelocity = 0;
        double magnormal = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            normalvelocity += normal[xyz]*Velocity<>::get<Lattice, Lattice::NDIM>(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq),xyz);
            magnormal += pow(normal[xyz],2);
        }

        magnormal = sqrt(magnormal);

        normalvelocity *= 1./magnormal;

        distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = (distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[idx]+normalvelocity*distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[idx])/(1+normalvelocity);
        //#pragma omp critical
        //{
        //std::cout<<k<<" "<<(distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[idx])<<std::endl;
        //}
        

    }    

}

template<class TTraits, class TDistributionType>
inline void Convective::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
    Lattice::communicateDistributionAllOld(distribution);

}
