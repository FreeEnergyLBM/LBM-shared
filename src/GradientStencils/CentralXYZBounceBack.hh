#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZBounceBack : GradientBase<Gradient, Cartesian> {
    template <class TTraits, class TParameter>
    inline double compute(int direction, int k);
};

template <class TTraits, class TParameter>
inline double CentralXYZBounceBack::compute(int direction, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    if (this->isBoundary<Lattice>(k)) return 0;

    double gradientsum = 0;
    const static auto& param = TParameter::template get<Lattice>();

    for (int idx = 1; idx < Stencil::Q; idx++) {
        if ((this->isBoundary<Lattice>(data.getNeighbor(k, idx)))) {
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (param[k]);

        } else {
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (param[data.getNeighbor(k, idx)]);
        }
    }
    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
}
