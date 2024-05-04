#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQ : GradientBase<Gradient, AllDirections> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k, int num = 0);
};

template <class TTraits, class TParameter>
inline double CentralQ::compute(const int direction, const int k, const int num) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) -
                  TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]],
                                                    num));
}