#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedQ : GradientBase<GradientMixed, AllDirections> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k, int num = 0);
};

template <class TTraits, class TParameter>
inline double MixedQ::compute(const int direction, const int k, int num) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    return 0.25 *
           (-TParameter::template get<Lattice>(
                data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction], num) +
            5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) -
            3 * TParameter::template get<Lattice>(k, num) -
            TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]],
                                              num));
}