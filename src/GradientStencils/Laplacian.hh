#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct LaplacianCentral : GradientBase<One> {

    template<class traits, class parameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class obj>
    using GradientType = Laplacian<obj,obj::instances>;
    
};


template<class traits, class parameter>
inline double LaplacianCentral::compute(const int direction, const int k, int num){
    
    using DataType = Data_Base<typename traits::Lattice, typename traits::Stencil>;

    double laplaciansum=0;

    for (int idx = 1; idx < traits::Stencil::Q; idx++) {

            laplaciansum +=  traits::Stencil::Weights[idx] * 2 * (parameter::template get<typename traits::Lattice>(DataType::getInstance().getNeighbors()[k * traits::Stencil::Q + idx],num) - parameter::template get<typename traits::Lattice>(k,num));

    }
    return 1.0 / (traits::Stencil::Cs2*traits::Lattice::m_DT*traits::Lattice::m_DT) * laplaciansum;
}
