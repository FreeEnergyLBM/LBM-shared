#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedXYZNoSolid : GradientBase<Cartesian> {

    template<class TTraits, class TParameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class TObj>
    using GradientType = GradientMixed<TObj,TObj::instances>;
    
};

template<class TTraits, class TParameter>
inline double MixedXYZNoSolid::compute(const int direction, const int k, int num){
        
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum=0;

    for (int idx = 1; idx < Stencil::Q; idx++) {

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])==1)) {

            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])!=1)) gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (2 *  TParameter::template get<Lattice>(k, num)
                                                                                       - 2 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

        }
        else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]] * Stencil::Q + Stencil::Opposites[direction]])==1)) {

            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)) gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q+  direction], num) 
                                                                                       - 4 * TParameter::template get<Lattice>(k, num));
            else gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q+  direction], num) 
                                                                                       - 3 * TParameter::template get<Lattice>(k, num)
                                                                                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

        }
        else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction]
                                                                                                                                                * Stencil::Q + direction], num)
                                                                                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                                                                                       - 4 * TParameter::template get<Lattice>(k, num));

        }
        else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])!=1)
                || (Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])!=1)) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction]
                                                                                                                                                * Stencil::Q+  direction], num)
                                                                                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                                                                                       - 3 * TParameter::template get<Lattice>(k, num)
                                                                                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

        }

    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
        
}