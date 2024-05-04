#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQBounceBack : GradientBase<Gradient, AllDirections> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k, int num = 0);
};

template <class TTraits, class TParameter>
inline double CentralQBounceBack::compute(const int direction, const int k, int num) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    if (this->isBoundary<Lattice>(k)) return 0;

    DataType& data = DataType::getInstance();
    const static auto& param = TParameter::template get<Lattice>();

    if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]))) {
        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) return 0;
        return 0.5 *
               (param[k * TParameter::instances + num] -
                param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]] * TParameter::instances +
                      num]);

    } else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
        return 0.5 * (param[data.getNeighbors()[k * Stencil::Q + direction] * TParameter::instances + num] -
                      param[k * TParameter::instances + num]);

    } else {
        return 0.5 *
               (param[data.getNeighbors()[k * Stencil::Q + direction] * TParameter::instances + num] -
                param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]] * TParameter::instances +
                      num]);
    }
}