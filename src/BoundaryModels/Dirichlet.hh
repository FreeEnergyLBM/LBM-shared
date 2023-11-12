#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class Dirichlet : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);

        inline void setInterfaceVal(double val) {mInterfaceVal=val;};

        inline void setInterfaceID(int id) {mInterfaceID[0]=id;};

        inline void setInterfaceID(const std::vector<int>& id) {mInterfaceID=id;};

    private:

        double mInterfaceVal;
        std::vector<int> mInterfaceID = {4};

};

template<class TTraits, class TDistributionType>
inline void Dirichlet::compute(TDistributionType& distribution, int k) {

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

            distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx)) + 2*TTraits::Stencil::Weights[idx]*mInterfaceVal;

        }    

}

template<class TTraits, class TDistributionType>
inline void Dirichlet::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}