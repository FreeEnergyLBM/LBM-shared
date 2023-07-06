#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZNoSolid : GradientBase<Cartesian> {
    
    template<class T_traits, class T_parameter>
    static inline double compute( int direction, int k, int num = 0);

    template<class T_obj>
    using GradientType = Gradient<T_obj,T_obj::instances>;
    
};

template<class T_traits, class T_parameter>
inline double CentralXYZNoSolid::compute(int direction, int k, int num) {
    
    using Lattice = typename T_traits::Lattice;
    using Stencil = typename T_traits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;

    for (int idx = 1; idx <Stencil::Q; idx++) {
        
        if ((Geometry<Lattice>::isSolid(data.getNeighbor(k,idx)))) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (T_parameter::template get<Lattice>(k, num));

        }
        else {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (T_parameter::template get<Lattice>(data.getNeighbor(k, idx), num));

        }
        
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;

}
