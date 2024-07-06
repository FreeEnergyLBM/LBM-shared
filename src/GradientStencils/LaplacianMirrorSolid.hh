#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct LaplacianCentralMirrorSolid : GradientBase<Laplacian, One> {  // FIX

    template <class TTraits, class TParameter>
    inline double compute(int direction, int k);
};

template <class TTraits, class TParameter>
inline double LaplacianCentralMirrorSolid::compute(const int direction, const int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (this->isBoundary<Lattice>(k)) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double laplaciansum = 0;
    const static auto& param = TParameter::template get<Lattice>();

    for (int idx = 1; idx < Stencil::Q; idx++) {
        if (!this->isBoundary<Lattice>(data.getNeighbor(k, idx))) {
            laplaciansum += Stencil::Weights[idx] * 2 * (param[data.getNeighbor(k, idx)] - param[k]);

        } else {
            const int& normalq =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(k, idx))
                              .NormalDirection)
                    ->second;

            double csolid = param[data.getNeighbor(data.getNeighbor(k, idx), normalq)];

            laplaciansum += Stencil::Weights[idx] * 2 * (csolid - param[k]);
        }
    }
    return 1.0 / (Stencil::Cs2 * Lattice::DT * Lattice::DT) * laplaciansum;
}
