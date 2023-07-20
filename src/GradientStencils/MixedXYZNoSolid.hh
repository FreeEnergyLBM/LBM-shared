#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedXYZNoSolid : GradientBase<Cartesian> {

    template<class T_traits, class T_parameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class T_obj>
    using GradientType = Gradient<T_obj,T_obj::instances>;
    
};

template<class T_traits, class T_parameter>
inline double MixedXYZNoSolid::compute(const int direction, const int k, int num){
        
    using Lattice = typename T_traits::Lattice;
    using Stencil = typename T_traits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum=0;

    for (int idx = 1; idx < Stencil::Q; idx++) {

        if ((Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + direction]))) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (2 *  T_parameter::template get<Lattice>(k, num)
                                                                                       - 2 * T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

        }
        else if ((Geometry<Lattice>::isSolid(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]]))) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (4 * T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q+  direction], num) 
                                                                                       - 3 * T_parameter::template get<Lattice>(k, num)
                                                                                       - T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

        }
        else if ((Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- T_parameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction]
                                                                                                                                                * Stencil::Q + direction], num)
                                                                                       + 5 * T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                                                                                       - 4 * T_parameter::template get<Lattice>(k, num));

        }
        else if ((!Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + direction]))
                || (!Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- T_parameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction]
                                                                                                                                                * Stencil::Q+  direction], num)
                                                                                       + 5 * T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                                                                                       - 3 * T_parameter::template get<Lattice>(k, num)
                                                                                       - T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

        }

    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
        
}