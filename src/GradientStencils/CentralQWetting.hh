#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQWetting : GradientBase<Gradient, AllDirections> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k);
};

template <class TTraits, class TParameter>
inline double CentralQWetting::compute(const int direction, const int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]))) {
        const int& normalq = TTraits::Stencil::QMap
                                 .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                           data.getNeighbor(k, direction))
                                           .NormalDirection)
                                 ->second;

        double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalq));

        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
            const int& normalqbackward =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(k, Stencil::Opposites[direction]))
                              .NormalDirection)
                    ->second;
            double csolidbackward =
                TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalqbackward));
            return 0.5 * ((csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))) -
                          (csolidbackward - 0.5 * this->mPrefactor * (csolidbackward - pow(csolidbackward, 2))));
        }
        return 0.5 *
               ((csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))) -
                TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]));

    } else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
        const int& normalq = TTraits::Stencil::QMap
                                 .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                           data.getNeighbor(k, Stencil::Opposites[direction]))
                                           .NormalDirection)
                                 ->second;

        double csolid = TParameter::template get<Lattice>(
            data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalq));

        return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]) -
                      (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))));

    } else {
        return 0.5 *
               (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]) -
                TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]));
    }
}
