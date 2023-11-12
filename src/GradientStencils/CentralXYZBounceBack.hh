#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZBounceBack : GradientBase<Cartesian> {
    
    template<class TTraits, class TParameter>
    static inline double compute( int direction, int k, int num = 0);

    template<class TObj>
    using GradientType = Gradient<TObj,TObj::instances>;
    
};

template<class TTraits, class TParameter>
inline double CentralXYZBounceBack::compute(int direction, int k, int num) {
    
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;

    for (int idx = 1; idx <Stencil::Q; idx++) {
        
        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k,idx))==1)) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(k, num));

        }
        else {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(data.getNeighbor(k, idx), num));

        }
        
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;

}
