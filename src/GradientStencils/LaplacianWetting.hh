#pragma once
#include <cmath>

#include "../Service.hh"
#include "GradientBase.hh"

struct LaplacianCentralWetting : GradientBase<Laplacian, One> {
    template <class TTraits, class TParameter>
    inline double compute(int direction, int k);
};

template <class TTraits, class TParameter>
inline double LaplacianCentralWetting::compute(int direction, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

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
            double factor = 0.5;  // Geometry<Lattice>::isCorner(data.getNeighbor(k, idx)) ? 1.0-sqrt(2)/2.0 : 0.5;
            laplaciansum += Stencil::Weights[idx] * 2 *
                            ((csolid - factor * this->mPrefactor * (csolid - pow(csolid, 2))) - param[k]);
        }
    }
    return 1.0 / (Stencil::Cs2 * Lattice::DT * Lattice::DT) * laplaciansum;
}
