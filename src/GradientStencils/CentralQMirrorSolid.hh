#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQMirrorSolid : GradientBase<Gradient, AllDirections> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k);
};

template <class TTraits, class TParameter>
inline double CentralQMirrorSolid::compute(const int direction, const int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();
    const static auto& param = TParameter::template get<Lattice>();

    if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]))) {
        const int& normalq = TTraits::Stencil::QMap
                                 .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                           data.getNeighbor(k, direction))
                                           .NormalDirection)
                                 ->second;

        double csolid = param[data.getNeighbor(data.getNeighbor(k, direction), normalq)];

        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
            const int& normalqbackward =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(k, Stencil::Opposites[direction]))
                              .NormalDirection)
                    ->second;
            double csolidbackward = param[data.getNeighbor(data.getNeighbor(k, direction), normalqbackward)];
            return 0.5 * (csolid - csolidbackward);
        }
        return 0.5 * (csolid - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]]);

    } else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
        const int& normalq = TTraits::Stencil::QMap
                                 .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                           data.getNeighbor(k, Stencil::Opposites[direction]))
                                           .NormalDirection)
                                 ->second;

        double csolid = param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalq)];

        return 0.5 * (param[data.getNeighbors()[k * Stencil::Q + direction]] - csolid);

    } else {
        return 0.5 * (param[data.getNeighbors()[k * Stencil::Q + direction]] -
                      param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]]);
    }
}
