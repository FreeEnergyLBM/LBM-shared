#pragma once
#include "BoundaryBase.hh"
#include "../Geometry.hh"
#include<iostream>


class FreeSlip : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        inline void setInterfaceID(int id) {mInterfaceID[0]=id;};

        inline void setInterfaceID(const std::vector<int>& id) {mInterfaceID=id;};

    private:

        double mInterfaceVal;
        std::vector<int> mInterfaceID = {1};

};

template<class TTraits, class TDistributionType>
inline void FreeSlip::compute(TDistributionType& distribution, int k) {

    for (int i : mInterfaceID){
        if(Geometry<typename TTraits::Lattice>::getBoundaryType(k) == i) goto runloop;
    }

    return;

    runloop:

        const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;
        const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection)->second;

        for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

            bool cont = true;

            for (int i : mInterfaceID){
                if(Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == i) goto dontapply;
            }

            cont = false;

            dontapply:
                if (cont) continue;

            std::array<int8_t,TTraits::Lattice::NDIM> newdir = {};

            newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[idx]-2*(int)normal[0]*(TTraits::Stencil::Ci_x[idx]==(int)normal[0]));
            if constexpr (TTraits::Lattice::NDIM>=2) newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[idx]-2*(int)normal[1]*(TTraits::Stencil::Ci_y[idx]==(int)normal[1]));
            if constexpr (TTraits::Lattice::NDIM>=3) newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[idx]-2*(int)normal[2]*(TTraits::Stencil::Ci_z[idx]==(int)normal[2]));

            const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;
            
            distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq),newidx))[newidx];//distribution.getPostCollisionDistribution(distribution.streamIndex(k, normalq),newidx);

        }    

}