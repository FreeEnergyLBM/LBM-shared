#pragma once
#include "BoundaryBase.hh"
#include "../Geometry.hh"
#include<iostream>


class BounceBack : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        inline void setInterfaceID(int id) {mInterfaceID[0]=id;};

        inline void setInterfaceID(const std::vector<int>& id) {mInterfaceID=id;};

    private:

        std::vector<int> mInterfaceID = {1};

};

template<class TTraits, class TDistributionType>
inline void BounceBack::compute(TDistributionType& distribution, int k) {

    for (int i : mInterfaceID){
        if(Geometry<typename TTraits::Lattice>::getBoundaryType(k) == i) goto runloop;
    }

    return;

    runloop:

        for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

            bool cont = true;

            for (int i : mInterfaceID){
                if(Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == i) goto dontapply;
            }

            cont = false;

            dontapply:
                if (cont) continue;
    
            distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx));

        }    

}