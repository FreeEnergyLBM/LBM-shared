#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZFourthOrder : GradientBase<Gradient, Cartesian> {
    template <class TTraits, class TParameter>
    inline double compute(int direction, int k);
};

template <class TTraits, class TParameter>
inline double CentralXYZFourthOrder::compute(int direction, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        const double& param1 = TParameter::template get<Lattice>(data.getNeighbor(k, idx));
        const double& param2 = TParameter::template get<Lattice>(data.getNeighbor(k, Stencil::Opposites[idx]));
        const double& param1_neighbor =
            TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx));
        const double& param2_neighbor = TParameter::template get<Lattice>(
            data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), Stencil::Opposites[idx]));

        gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] *
                       (-param1_neighbor + 8 * param1 - 8 * param2 + param2_neighbor);
    }

    return 1.0 / (12.0 * Stencil::Cs2 * Lattice::DT) * gradientsum;
}
