#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedQBounceBack : GradientBase<GradientMixed, AllDirections> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k, int num = 0);
};

template <class TTraits, class TParameter>
inline double MixedQBounceBack::compute(const int direction, const int k, int num) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    if (this->isBoundary<Lattice>(k)) return 0;

    DataType& data = DataType::getInstance();

    const static auto& param = TParameter::template get<Lattice>();

    if ((!this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction])) &&
        (!this->isBoundary<Lattice>(
            data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction])) &&
        (!this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
        return 0.25 *
               (-param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction] *
                           TParameter::instances +
                       num] +
                5 * param[data.getNeighbors()[k * Stencil::Q + direction] * TParameter::instances + num] -
                3 * param[k * TParameter::instances + num] -
                param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]] * TParameter::instances +
                      num]);

    } else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]))) {
        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
            return 0;
        }

        return 0.25 *
               (2 * param[k * TParameter::instances + num] -
                2 * param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]] * TParameter::instances +
                          num]);

    } else if ((this->isBoundary<Lattice>(
                   data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction]))) {
        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
            return 0.25 * (4 * param[data.getNeighbors()[k * Stencil::Q + direction] * TParameter::instances + num] -
                           4 * param[k * TParameter::instances + num]);
        }

        return 0.25 *
               (4 * param[data.getNeighbors()[k * Stencil::Q + direction] * TParameter::instances + num] -
                3 * param[k * TParameter::instances + num] -
                param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]] * TParameter::instances +
                      num]);

    }

    else {
        return 0.25 *
               (-param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction] *
                           TParameter::instances +
                       num] +
                5 * param[data.getNeighbors()[k * Stencil::Q + direction] * TParameter::instances + num] -
                4 * param[k * TParameter::instances + num]);
    }
}