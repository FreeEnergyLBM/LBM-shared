#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZFourthOrder : GradientBase<Cartesian> {

    template<class T_traits, class T_parameter>
    static inline double compute(int direction, int k, int num = 0);

    template<class T_obj>
    using GradientType = Gradient<T_obj,T_obj::instances>;
    
};

template<class T_traits, class T_parameter>
inline double CentralXYZFourthOrder::compute(int direction, int k, int num){
    
    using Lattice = typename T_traits::Lattice;
    using Stencil = typename T_traits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;

    for (int idx = 1; idx <Stencil::Q; idx++) {

            const double& param1 = T_parameter::template get<Lattice>(data.getNeighbor(k, idx), num);
            const double& param2 = T_parameter::template get<Lattice>(data.getNeighbor(k, Stencil::Opposites[idx]), num);
            const double& param1_neighbor = T_parameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx), num);
            const double& param2_neighbor = T_parameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), Stencil::Opposites[idx]), num);

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (- param1_neighbor + 8 * param1 - 8 * param2 + param2_neighbor);
        
    }

    return 1.0 / (12.0 * Stencil::Cs2 * Lattice::DT) * gradientsum;

}
