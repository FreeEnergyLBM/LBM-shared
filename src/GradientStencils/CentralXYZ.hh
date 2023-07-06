#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZ : GradientBase<Cartesian> {

    template<class traits, class parameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class obj>
    using GradientType = Gradient<obj,obj::instances>;
    
};

template<class traits, class parameter>
inline double CentralXYZ::compute(const int direction, const int k, int num){
    
    using DataType = Data_Base<typename traits::Lattice, typename traits::Stencil>;

    double gradientsum=0;

    for (int idx = 1; idx <traits::Stencil::Q; idx++) {

            gradientsum += 0.5*traits::Stencil::Weights[idx] * traits::Stencil::Ci_xyz(direction)[idx] *(parameter::template get<typename traits::Lattice>(DataType::getInstance().getNeighbors()[k * traits::Stencil::Q+idx],num)-parameter::template get<typename traits::Lattice>(DataType::getInstance().getNeighbors()[k * traits::Stencil::Q + traits::Stencil::Opposites[idx]],num)); //TEMPORARY
        
    }

    return 1.0 / (traits::Stencil::Cs2*traits::Lattice::m_DT) * gradientsum;

}
