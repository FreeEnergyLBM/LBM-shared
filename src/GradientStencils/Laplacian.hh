#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct LaplacianCentral : GradientBase<One> {

    template<class TTraits, class TParameter>
    inline double compute(int direction, int k, int num = 0);

    template<class TObj>
    using GradientType = Laplacian<TObj,TObj::instances>;
    
};


template<class TTraits, class TParameter>
inline double LaplacianCentral::compute(int direction, int k, int num){

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double laplaciansum=0;

    for (int idx = 1; idx < Stencil::Q; idx++) {

            laplaciansum +=  Stencil::Weights[idx] * 2 * (TParameter::template get<Lattice>(data.getNeighbor(k, idx), num) - TParameter::template get<Lattice>(k, num));

    }
    return 1.0 / (Stencil::Cs2 * Lattice::DT * Lattice::DT) * laplaciansum;
}
