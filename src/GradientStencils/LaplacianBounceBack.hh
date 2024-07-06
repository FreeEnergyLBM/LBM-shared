#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct LaplacianCentralBounceBack : GradientBase<Laplacian, One> {  // FIX

    template <class TTraits, class TParameter>
    inline double compute(int direction, int k);
};

template <class TTraits, class TParameter>
inline double LaplacianCentralBounceBack::compute(const int direction, const int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double laplaciansum = 0;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        if (!this->isBoundary<Lattice>(data.getNeighbor(k, idx))) {
            laplaciansum +=
                Stencil::Weights[idx] * 2 *
                (TParameter::template get<Lattice>(data.getNeighbor(k, idx)) - TParameter::template get<Lattice>(k));
        }
    }
    return 1.0 / (Stencil::Cs2 * Lattice::DT * Lattice::DT) * laplaciansum;
}
