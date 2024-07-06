#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedQWetting : GradientBase<GradientMixed, AllDirections> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k);
};

template <class TTraits, class TParameter>
inline double MixedQWetting::compute(const int direction, const int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (this->isBoundary<Lattice>(k)) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    if ((this->isBoundary<Lattice>(data.getNeighbor(k, direction)))) {
        if ((this->isBoundary<Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])))) {
            const std::array<int8_t, TTraits::Lattice::NDIM>& normal =
                BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                    data.getNeighbor(k, direction))
                    .NormalDirection;
            const int& normalqforward =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(k, direction))
                              .NormalDirection)
                    ->second;
            const int& normalqbackward =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(k, Stencil::Opposites[direction]))
                              .NormalDirection)
                    ->second;

            std::array<int8_t, TTraits::Lattice::NDIM> newdir = {};

            newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[direction] +
                                 2 * (int)normal[0] * (TTraits::Stencil::Ci_x[direction] == -(int)normal[0]));
            if constexpr (TTraits::Lattice::NDIM >= 2)
                newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[direction] +
                                     2 * (int)normal[1] * (TTraits::Stencil::Ci_y[direction] == -(int)normal[1]));
            if constexpr (TTraits::Lattice::NDIM >= 3)
                newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[direction] +
                                     2 * (int)normal[2] * (TTraits::Stencil::Ci_z[direction] == -(int)normal[2]));

            const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;

            double csolid =
                TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalqforward));
            double csolid2 = TParameter::template get<Lattice>(data.getNeighbor(
                data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalqforward]), newidx), newidx));
            double csolid3 = TParameter::template get<Lattice>(
                data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalqbackward));
            return 0.25 * (-((csolid2 - 0.5 * this->mPrefactor * (csolid2 - pow(csolid2, 2)))) +
                           5 * ((csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2)))) -
                           3 * TParameter::template get<Lattice>(k) -
                           (csolid3 - 0.5 * this->mPrefactor * (csolid3 - pow(csolid3, 2))));
        } else {
            const std::array<int8_t, TTraits::Lattice::NDIM>& normal =
                BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                    data.getNeighbor(k, direction))
                    .NormalDirection;
            const int& normalq =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(k, direction))
                              .NormalDirection)
                    ->second;

            std::array<int8_t, TTraits::Lattice::NDIM> newdir = {};

            newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[direction] +
                                 2 * (int)normal[0] * (TTraits::Stencil::Ci_x[direction] == -(int)normal[0]));
            if constexpr (TTraits::Lattice::NDIM >= 2)
                newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[direction] +
                                     2 * (int)normal[1] * (TTraits::Stencil::Ci_y[direction] == -(int)normal[1]));
            if constexpr (TTraits::Lattice::NDIM >= 3)
                newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[direction] +
                                     2 * (int)normal[2] * (TTraits::Stencil::Ci_z[direction] == -(int)normal[2]));

            const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;

            double csolid =
                TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalq));
            double csolid2 = TParameter::template get<Lattice>(
                data.getNeighbor(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalq]), newidx), newidx));

            return 0.25 * (-((csolid2 - 0.5 * this->mPrefactor * (csolid2 - pow(csolid2, 2)))) +
                           5 * ((csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2)))) -
                           3 * TParameter::template get<Lattice>(k) -
                           TParameter::template get<Lattice>(
                               data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]));
        }

    } else if (this->isBoundary<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), direction))) {
        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
            const int& normalqforward =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(data.getNeighbor(k, direction), direction))
                              .NormalDirection)
                    ->second;
            double csolidforward = TParameter::template get<Lattice>(
                data.getNeighbor(data.getNeighbor(data.getNeighbor(k, direction), direction), normalqforward));
            const int& normalqbackward =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(k, Stencil::Opposites[direction]))
                              .NormalDirection)
                    ->second;

            double csolidbackward = TParameter::template get<Lattice>(
                data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalqbackward));

            return 0.25 * (-(csolidforward - 0.5 * this->mPrefactor * (csolidforward - pow(csolidforward, 2))) +
                           5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]) -
                           3 * TParameter::template get<Lattice>(k) -
                           (csolidbackward - 0.5 * this->mPrefactor * (csolidbackward - pow(csolidbackward, 2))));
        }

        const int& normalq = TTraits::Stencil::QMap
                                 .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                           data.getNeighbor(data.getNeighbor(k, direction), direction))
                                           .NormalDirection)
                                 ->second;
        double csolid = TParameter::template get<Lattice>(
            data.getNeighbor(data.getNeighbor(data.getNeighbor(k, direction), direction), normalq));

        return 0.25 *
               (-(csolid) + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]) -
                3 * TParameter::template get<Lattice>(k) -
                TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]));

    }

    else if ((this->isBoundary<Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])) == 1)) {
        const int& normalq = TTraits::Stencil::QMap
                                 .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                           data.getNeighbor(k, Stencil::Opposites[direction]))
                                           .NormalDirection)
                                 ->second;
        double csolid = TParameter::template get<Lattice>(
            data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalq));

        return 0.25 *
               (-TParameter::template get<Lattice>(
                    data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction]) +
                5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]) -
                3 * TParameter::template get<Lattice>(k) -
                ((csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2)))));

    } else {
        return 0.25 *
               (-TParameter::template get<Lattice>(
                    data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction]) +
                5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]) -
                3 * TParameter::template get<Lattice>(k) -
                TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]));
    }
}
