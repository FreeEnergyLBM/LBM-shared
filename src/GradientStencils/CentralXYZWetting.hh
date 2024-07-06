#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZWetting : GradientBase<Gradient, Cartesian> {
    template <class TTraits, class TParameter>
    inline double compute(int direction, int k);
};

template <class TTraits, class TParameter>
inline double CentralXYZWetting::compute(int direction, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        if ((this->isBoundary<Lattice>(data.getNeighbor(k, idx)))) {
            const int& normalq =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(k, idx))
                              .NormalDirection)
                    ->second;

            double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), normalq));

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] *
                           (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2)));

        } else {
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] *
                           (TParameter::template get<Lattice>(data.getNeighbor(k, idx)));
        }
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
}
