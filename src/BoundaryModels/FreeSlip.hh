#pragma once
#include "BoundaryBase.hh"
#include "../Geometry.hh"
#include<iostream>


class FreeSlip : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

    private:

};

template<class TTraits, class TDistributionType>
inline void FreeSlip::compute(TDistributionType& distribution, int k) {

    if (Geometry<typename TTraits::Lattice>::getBoundaryType(k) != 1) return;

    std::vector<int8_t>& normal = BoundaryLabels<>::get<typename TTraits::Lattice>(k).NormalDirection;
    const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<>::get<typename TTraits::Lattice>(k).NormalDirection)->second;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

        if(Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == 0 ) {

            std::vector<int8_t> newdir(TTraits::Lattice::NDIM,0);

            newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[idx]-2*(int)normal[0]*(TTraits::Stencil::Ci_x[idx]!=0));
            if constexpr (TTraits::Lattice::NDIM>=2) newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[idx]-2*(int)normal[1]*(TTraits::Stencil::Ci_y[idx]!=0));
            if constexpr (TTraits::Lattice::NDIM>=3) newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[idx]-2*(int)normal[2]*(TTraits::Stencil::Ci_z[idx]!=0));

            const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;
            
            distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = distribution.getDistributionPointer(k)[newidx];
        
        }

    }    

}