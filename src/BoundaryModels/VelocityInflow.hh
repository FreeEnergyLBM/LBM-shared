#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class VelocityInflow : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);

        inline void setWallVelocity(const std::vector<double>& momentum) {mWallMomentum=momentum;};

        inline void setInterfaceID(int id) {mInterfaceID[0]=id;};

        inline void setInterfaceID(const std::vector<int>& id) {mInterfaceID=id;};

    private:

        std::vector<double> mWallMomentum;
        std::vector<int> mInterfaceID = {3};

};

template<class TTraits, class TDistributionType>
inline void VelocityInflow::compute(TDistributionType& distribution, int k) { //CHANGE THIS SO YOU DONT NEED TO COMMUNICATE

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

            double cidotmomentum = 0;
            for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
                cidotmomentum += TTraits::Stencil::Ci_xyz(xyz)[idx] * mWallMomentum[xyz];
            }
            distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx)) - 2*TTraits::Stencil::Weights[idx]*cidotmomentum/TTraits::Stencil::Cs2;

        }    

}

template<class TTraits, class TDistributionType>
inline void VelocityInflow::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}