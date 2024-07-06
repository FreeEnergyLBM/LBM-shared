#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedQBounceBack : GradientBase<GradientMixed, AllDirections> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k);
};

template <class TTraits, class TParameter>
inline double MixedQBounceBack::compute(const int direction, const int k) {
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
               (-param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction]] +
                5 * param[data.getNeighbors()[k * Stencil::Q + direction]] - 3 * param[k] -
                param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]]);

    } else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]))) {
        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
            return 0;
        }

        return 0.25 * (2 * param[k] - 2 * param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]]);

    } else if ((this->isBoundary<Lattice>(
                   data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction]))) {
        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
            return 0.25 * (4 * param[data.getNeighbors()[k * Stencil::Q + direction]] - 4 * param[k]);
        }

        return 0.25 * (4 * param[data.getNeighbors()[k * Stencil::Q + direction]] - 3 * param[k] -
                       param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]]);

    }

    else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
        return 0.25 *
               (-param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction]] +
                5 * param[data.getNeighbors()[k * Stencil::Q + direction]] - 4 * param[k]);
    } else
        return 0;
}
