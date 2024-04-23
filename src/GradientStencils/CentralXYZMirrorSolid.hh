#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZMirrorSolid : GradientBase<Gradient,Cartesian> {
    
    template<class TTraits, class TParameter>
    inline double compute( int direction, int k, int num = 0);
    
};

template<class TTraits, class TParameter>
inline double CentralXYZMirrorSolid::compute(int direction, int k, int num) {
    
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    //if (this->isBoundary<Lattice>(k) == 4) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;
    const static auto& param = TParameter::template get<Lattice>();
    const static auto& boundary = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>();
    for (int idx = 1; idx <Stencil::Q; idx++) {
        
        if ((this->isBoundary<Lattice>(data.getNeighbor(k,idx)))) {
            
                const int& normalq = TTraits::Stencil::QMap.find(boundary[data.getNeighbor(k, idx)].NormalDirection)->second;

                double csolid = param[data.getNeighbor(data.getNeighbor(k, idx), normalq)*TParameter::instances + num];

                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * csolid;

        }
        else {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (param[data.getNeighbor(k, idx)*TParameter::instances + num]);

        }
        
    }
    
    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;

}
