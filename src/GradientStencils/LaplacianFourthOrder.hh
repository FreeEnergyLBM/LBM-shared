#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct LaplacianCentralFourthOrder : GradientBase<Laplacian,One> {

    template<class TTraits, class TParameter>
    inline double compute(int direction, int k, int num = 0);
    
};


template<class TTraits, class TParameter>
inline double LaplacianCentralFourthOrder::compute(int direction, int k, int num){

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double laplaciansum=0;

    for (int idx = 1; idx < Stencil::Q; idx++) {

            laplaciansum +=  Stencil::Weights[idx] * 2 * (- TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx), num) + 16 * TParameter::template get<Lattice>(data.getNeighbor(k, idx), num) - 15 * TParameter::template get<Lattice>(k, num));

    }
    return 1.0 / (12.0 * Stencil::Cs2 * Lattice::DT * Lattice::DT) * laplaciansum;
}
