#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class ExtrapolationOutflow : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);

        inline void setInterfaceID(int id) {mInterfaceID[0]=id;};

        inline void setInterfaceID(const std::vector<int>& id) {mInterfaceID=id;};

    private:

        double mInterfaceVal;
        std::vector<int> mInterfaceID = {4};

};

template<class TTraits, class TDistributionType>
inline void ExtrapolationOutflow::compute(TDistributionType& distribution, int k) { //CHANGE THIS SO YOU DONT NEED TO COMMUNICATE
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    for (int i : mInterfaceID){
        if(Geometry<typename TTraits::Lattice>::getBoundaryType(k) == i) goto runloop;
    }

    return;

    runloop:

        const int& normalq = Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;

        for (int idx = 1; idx < Stencil::Q; idx++) {
            
            bool cont = true;

            for (int i : mInterfaceID){
                if(Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == i) goto dontapply;
            }

            cont = false;

            dontapply:
                if (cont) continue;

            double cidotnormal = 0;
            for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
                cidotnormal += TTraits::Stencil::Ci_xyz(xyz)[idx]*BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection[xyz];
            }
            
            if(cidotnormal!=0) {
                
                distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = (4.0*distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[idx] - distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq), normalq))[idx])/3.0;

            
            }
            else{
                distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = distribution.getEquilibriumPointer(distribution.streamIndex(k, normalq))[idx];//TTraits::Stencil::Weights
            }

        }    

}

template<class TTraits, class TDistributionType>
inline void ExtrapolationOutflow::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
    Lattice::communicateDistributionAllEquilibrium(distribution);

}
