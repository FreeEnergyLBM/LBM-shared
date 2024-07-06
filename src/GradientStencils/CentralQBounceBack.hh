#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQBounceBack : GradientBase<Gradient, AllDirections> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k);
};

template <class TTraits, class TParameter>
inline double CentralQBounceBack::compute(const int direction, const int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    if (this->isBoundary<Lattice>(k)) return 0;

    DataType& data = DataType::getInstance();
    const static auto& param = TParameter::template get<Lattice>();

    if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]))) {
        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) return 0;
        return 0.5 * (param[k] - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]]);

    } else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
        return 0.5 * (param[data.getNeighbors()[k * Stencil::Q + direction]] - param[k]);

    } else {
        return 0.5 * (param[data.getNeighbors()[k * Stencil::Q + direction]] -
                      param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]]);
    }
}
