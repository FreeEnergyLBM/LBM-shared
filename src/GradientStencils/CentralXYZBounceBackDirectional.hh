#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZBounceBackDirectional : GradientBase<Gradient, Cartesian> {
    template <class TTraits, class TParameter, int TNdir>
    inline double compute(int direction, int k, int dir);
};

template <class TTraits, class TParameter, int TNdir>
inline double CentralXYZBounceBackDirectional::compute(int graddir, int k, int dir) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    if (this->isBoundary<Lattice>(k)) return 0;

    double gradientsum = 0;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        if ((this->isBoundary<Lattice>(data.getNeighbor(k, idx)))) {
            gradientsum +=
                Stencil::Weights[idx] * Stencil::Ci_xyz(graddir)[idx] * (-TParameter::template get<Lattice,TNdir>(data.getNeighbor(k, idx), dir));

        } else {
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(graddir)[idx] *
                           (TParameter::template get<Lattice,TNdir>(data.getNeighbor(k, idx), dir));
        }
    }
    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
}
