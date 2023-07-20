#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedXYZ : GradientBase<Cartesian> {

    template<class T_traits, class T_parameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class T_obj>
    using GradientType = Gradient<T_obj,T_obj::instances>;
    
};

template<class T_traits, class T_parameter>
inline double MixedXYZ::compute(const int direction, const int k, int num){
    
    using Lattice = typename T_traits::Lattice;
    using Stencil = typename T_traits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum=0;

    for (int idx = 1; idx < Stencil::Q; idx++) {

        gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (- T_parameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx], num) 
                                                                                  + 5 * T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num) 
                                                                                  - 3 * T_parameter::template get<Lattice>(k, num)
                                                                                  - T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));
        
    }

    return 0.25 / (Stencil::Cs2*Lattice::DT) * gradientsum;

}