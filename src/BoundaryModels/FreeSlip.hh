#pragma once
#include "BoundaryBase.hh"
#include<iostream>


class FreeSlip : public BoundaryBase {
    public:

        FreeSlip() { this->setInterfaceID(1); }

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);

        template<class TTraits>
        inline void communicate(){};

};

template<class TTraits, class TDistributionType>
inline void FreeSlip::compute(TDistributionType& distribution, int k) {

    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;

    const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;
    int normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection)->second;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

        if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;

        if (Geometry<typename TTraits::Lattice>::isCorner(k)) {

            std::array<int8_t,TTraits::Lattice::NDIM> cinorm = {};

            cinorm[0] = (int8_t)((int)normal[0]*(TTraits::Stencil::Ci_x[idx]==(int)normal[0]));
            if constexpr (TTraits::Lattice::NDIM>=2) cinorm[1] = (int8_t)((int)normal[1]*(TTraits::Stencil::Ci_y[idx]==(int)normal[1]));
            if constexpr (TTraits::Lattice::NDIM>=3) cinorm[2] = (int8_t)((int)normal[2]*(TTraits::Stencil::Ci_z[idx]==(int)normal[2]));

            normalq = TTraits::Stencil::QMap.find(cinorm)->second;

        }

        std::array<int8_t,TTraits::Lattice::NDIM> newdir = {};

        newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[idx]-2*(int)normal[0]*(TTraits::Stencil::Ci_x[idx]==(int)normal[0]));
        if constexpr (TTraits::Lattice::NDIM>=2) newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[idx]-2*(int)normal[1]*(TTraits::Stencil::Ci_y[idx]==(int)normal[1]));
        if constexpr (TTraits::Lattice::NDIM>=3) newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[idx]-2*(int)normal[2]*(TTraits::Stencil::Ci_z[idx]==(int)normal[2]));

        const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;
        
        distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq),newidx))[newidx];//distribution.getPostCollisionDistribution(distribution.streamIndex(k, normalq),newidx);

    }    

}

template<class TTraits, class TDistributionType>
inline void FreeSlip::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}
