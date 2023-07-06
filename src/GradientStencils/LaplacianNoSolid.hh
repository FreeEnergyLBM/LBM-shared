#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct LaplacianCentralNoSolid : GradientBase<One> { //FIX

    template<class T_traits, class T_parameter>
    static inline double compute(int direction, int k, int num = 0);

    template<class T_obj>
    using GradientType = Laplacian<T_obj,T_obj::instances>;
    
};

template<class T_traits, class T_parameter>
inline double LaplacianCentralNoSolid::compute(const int direction, const int k, int num){
   
    using Lattice = typename T_traits::Lattice;
    using Stencil = typename T_traits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double laplaciansum=0;

    for (int idx = 1; idx <Stencil::Q; idx++) {

            if(!Geometry<Lattice>::isSolid(data.getNeighbor(k, idx))) {

                laplaciansum +=  Stencil::Weights[idx] * 2 * (T_parameter::template get<Lattice>(data.getNeighbor(k, idx), num) - T_parameter::template get<Lattice>(k, num));

            }

    }
    return 1.0 / (Stencil::Cs2 * Lattice::DT * Lattice::DT) * laplaciansum;
}
