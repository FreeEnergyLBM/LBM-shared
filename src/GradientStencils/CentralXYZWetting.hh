#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZWetting : WettingGradient<Cartesian> {

    template<class TTraits, class TParameter>
    inline double compute(int direction, int k, int num = 0);

    template<class TObj>
    using GradientType = Gradient<TObj,TObj::instances>;
    
};
/*
template<class TTraits, class TParameter>
inline double CentralXYZWetting::compute(int direction, int k, int num){
    
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    //if (this->isBoundary<Lattice>(k) == 4) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;

    for (int idx = 1; idx <Stencil::Q; idx++) {
        
        if ((this->isBoundary<Lattice>(data.getNeighbor(k,idx))==1)) {

            double csolid = TParameter::template get<Lattice>(k, num);

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2)));

        }
        else {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(data.getNeighbor(k, idx), num));

        }
        
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;

}
*/


template<class TTraits, class TParameter>
inline double CentralXYZWetting::compute(int direction, int k, int num) {
    
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    //if (this->isBoundary<Lattice>(k) == 4) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;

    for (int idx = 1; idx <Stencil::Q; idx++) {
        
        if ((this->isBoundary<Lattice>(data.getNeighbor(k,idx)))) {
            
                const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;

                double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), normalq), num);

                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2)));

        }
        else {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(data.getNeighbor(k, idx), num));

        }
        
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;

}
