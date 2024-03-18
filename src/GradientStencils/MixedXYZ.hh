#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedXYZ : GradientBase<Cartesian> {

    template<class TTraits, class TParameter>
    inline double compute(const int direction, const int k, int num = 0);

    template<class TObj>
    using GradientType = GradientMixed<TObj,TObj::instances>;
    
};

template<class TTraits, class TParameter>
inline double MixedXYZ::compute(const int direction, const int k, int num){
    
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum=0;

    for (int idx = 1; idx < Stencil::Q; idx++) {

        gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx], num) 
                                                                                  + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num) 
                                                                                  - 3 * TParameter::template get<Lattice>(k, num)
                                                                                  - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));
        
    }

    return 0.25 / (Stencil::Cs2*Lattice::DT) * gradientsum;

}