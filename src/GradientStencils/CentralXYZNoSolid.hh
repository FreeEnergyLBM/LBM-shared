#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZNoSolid : GradientBase<Cartesian> {
    
    template<class traits, class parameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class obj>
    using GradientType = Gradient<obj,obj::instances>;
    
};

template<class traits, class parameter>
inline double CentralXYZNoSolid::compute(const int direction, const int k, int num) {
    
    using DataType = Data_Base<typename traits::Lattice, typename traits::Stencil>;

    double gradientsum=0;

    for (int idx = 1; idx <traits::Stencil::Q; idx++) {
        
        if ((Geometry<typename traits::Lattice>::isSolid(DataType::getInstance().getNeighbors()[k * traits::Stencil::Q + idx]))) {

            gradientsum += traits::Stencil::Weights[idx] * traits::Stencil::Ci_xyz(direction)[idx] * (parameter::template get<typename traits::Lattice>(k,num));

        }
        else {

            gradientsum += traits::Stencil::Weights[idx] * traits::Stencil::Ci_xyz(direction)[idx] * (parameter::template get<typename traits::Lattice>(DataType::getInstance().getNeighbors()[k * traits::Stencil::Q+idx],num));

        }
        
    }

    return 1.0 / (traits::Stencil::Cs2*traits::Lattice::m_DT) * gradientsum;

}
