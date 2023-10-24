#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class VelocityInflow : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicatePrecompute(){};

        template<class TTraits, class TDistributionType>
        inline void communicatePrecompute(TDistributionType& mDistribution);

        inline void setWallVelocity(const std::vector<double>& momentum) {mWallMomentum=momentum;};

        inline void setInterfaceID(int id) {mInterfaceID=id;};

    private:

        std::vector<double> mWallMomentum;
        int mInterfaceID;

};

template<class TTraits, class TDistributionType>
inline void VelocityInflow::compute(TDistributionType& distribution, int k) { //CHANGE THIS SO YOU DONT NEED TO COMMUNICATE

    if (Geometry<typename TTraits::Lattice>::getBoundaryType(k) != mInterfaceID) return;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

        if(Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == 0 ) {
            
            //distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -distribution.getDistributionOldPointer(distribution.streamIndex(k, idx))[distribution.getOpposite(idx)] + 2*TTraits::Stencil::Weights[idx]*mInterfaceVal;
            double cidotmomentum = 0;
            for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
                cidotmomentum += TTraits::Stencil::Ci_xyz(xyz)[idx] * mWallMomentum[xyz];
            }
            distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = distribution.getDistributionPointer(k)[distribution.getOpposite(idx)] - 2*TTraits::Stencil::Weights[idx]*cidotmomentum/TTraits::Stencil::Cs2;//TTraits::Stencil::Weights[idx]*
            //distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -distribution.getDistributionPointer(k)[distribution.getOpposite(idx)] + 2*TTraits::Stencil::Weights[idx]*mInterfaceVal;
        
        }

    }    

}

template<class TTraits, class TDistributionType>
inline void VelocityInflow::communicatePrecompute(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}