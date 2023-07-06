#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZ : GradientBase<Cartesian> {

    template<class T_traits, class T_parameter>
    static inline double compute(int direction, int k, int num = 0);

    template<class T_obj>
    using GradientType = Gradient<T_obj,T_obj::instances>;
    
};

template<class T_traits, class T_parameter>
inline double CentralXYZ::compute(int direction, int k, int num){
    
    using Lattice = typename T_traits::Lattice;
    using Stencil = typename T_traits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;

    for (int idx = 1; idx <Stencil::Q; idx++) {

            const double& param1 = T_parameter::template get<Lattice>(data.getNeighbor(k, idx), num);
            const double& param2 = T_parameter::template get<Lattice>(data.getNeighbor(k, Stencil::Opposites[idx]), num);

            gradientsum += 0.5 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (param1 - param2);
        
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;

}
