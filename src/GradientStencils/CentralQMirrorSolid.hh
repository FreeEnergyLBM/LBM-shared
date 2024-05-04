#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQMirrorSolid : GradientBase<Gradient, AllDirections> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k, int num = 0);
};

template <class TTraits, class TParameter>
inline double CentralQMirrorSolid::compute(const int direction, const int k, int num) {
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

        double csolid = param[data.getNeighbor(data.getNeighbor(k, direction), normalq) * TParameter::instances + num];

        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
            const int& normalqbackward =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(k, Stencil::Opposites[direction]))
                              .NormalDirection)
                    ->second;
            double csolidbackward =
                param[data.getNeighbor(data.getNeighbor(k, direction), normalqbackward) * TParameter::instances + num];
            return 0.5 * (csolid - csolidbackward);
        }
        return 0.5 *
               (csolid -
                param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]] * TParameter::instances +
                      num]);

    } else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
        const int& normalq = TTraits::Stencil::QMap
                                 .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                           data.getNeighbor(k, Stencil::Opposites[direction]))
                                           .NormalDirection)
                                 ->second;

        double csolid = param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalq) *
                                  TParameter::instances +
                              num];

        return 0.5 * (param[data.getNeighbors()[k * Stencil::Q + direction] * TParameter::instances + num] - csolid);

    } else {
        return 0.5 *
               (param[data.getNeighbors()[k * Stencil::Q + direction] * TParameter::instances + num] -
                param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]] * TParameter::instances +
                      num]);
    }
}