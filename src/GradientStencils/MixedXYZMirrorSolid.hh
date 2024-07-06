#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedXYZMirrorSolid : GradientBase<GradientMixed, Cartesian> {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k);
};

template <class TTraits, class TParameter>
inline double MixedXYZMirrorSolid::compute(const int direction, const int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    static DataType& data = DataType::getInstance();
    const static auto& param = TParameter::template get<Lattice>();
    double gradientsum = 0;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        if ((!this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + idx])) &&
            (!this->isBoundary<Lattice>(
                data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx])) &&
            (!this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 *
                           (-param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx]] +
                            5 * param[data.getNeighbors()[k * Stencil::Q + idx]] - 3 * param[k] -
                            param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]]);

        } else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + idx]))) {
            if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {
                const std::array<int8_t, TTraits::Lattice::NDIM>& normal =
                    BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                        data.getNeighbor(k, idx))
                        .NormalDirection;
                const int& normalq =
                    TTraits::Stencil::QMap
                        .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                  data.getNeighbor(k, idx))
                                  .NormalDirection)
                        ->second;
                const int& normalqbackward =
                    TTraits::Stencil::QMap
                        .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                  data.getNeighbor(k, Stencil::Opposites[idx]))
                                  .NormalDirection)
                        ->second;

                std::array<int8_t, TTraits::Lattice::NDIM> newdir = {};

                newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[idx] +
                                     2 * (int)normal[0] * (TTraits::Stencil::Ci_x[idx] == -(int)normal[0]));
                if constexpr (TTraits::Lattice::NDIM >= 2)
                    newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[idx] +
                                         2 * (int)normal[1] * (TTraits::Stencil::Ci_y[idx] == -(int)normal[1]));
                if constexpr (TTraits::Lattice::NDIM >= 3)
                    newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[idx] +
                                         2 * (int)normal[2] * (TTraits::Stencil::Ci_z[idx] == -(int)normal[2]));

                const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;

                double csolid = param[data.getNeighbor(data.getNeighbor(k, idx), normalq)];
                double csolid2 = param[data.getNeighbor(
                    data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalq]), newidx), newidx)];
                double csolid3 = param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), normalqbackward)];

                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 *
                               (-(csolid2) + 5 * (csolid)-3 * param[k] - csolid3);

            } else {
                const std::array<int8_t, TTraits::Lattice::NDIM>& normal =
                    BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                        data.getNeighbor(k, idx))
                        .NormalDirection;
                const int& normalq =
                    TTraits::Stencil::QMap
                        .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                  data.getNeighbor(k, idx))
                                  .NormalDirection)
                        ->second;

                std::array<int8_t, TTraits::Lattice::NDIM> newdir = {};

                newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[idx] +
                                     2 * (int)normal[0] * (TTraits::Stencil::Ci_x[idx] == -(int)normal[0]));
                if constexpr (TTraits::Lattice::NDIM >= 2)
                    newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[idx] +
                                         2 * (int)normal[1] * (TTraits::Stencil::Ci_y[idx] == -(int)normal[1]));
                if constexpr (TTraits::Lattice::NDIM >= 3)
                    newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[idx] +
                                         2 * (int)normal[2] * (TTraits::Stencil::Ci_z[idx] == -(int)normal[2]));

                const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;

                double csolid = param[data.getNeighbor(data.getNeighbor(k, idx), normalq)];
                double csolid2 = param[data.getNeighbor(
                    data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalq]), newidx), newidx)];

                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 *
                               (-(csolid2) + 5 * (csolid)-3 * param[k] -
                                param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]]);
            }

        } else if ((this->isBoundary<Lattice>(
                       data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx]))) {
            if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {
                const int& normalqforward =
                    TTraits::Stencil::QMap
                        .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                  data.getNeighbor(data.getNeighbor(k, idx), idx))
                                  .NormalDirection)
                        ->second;
                double csolidforward =
                    param[data.getNeighbor(data.getNeighbor(data.getNeighbor(k, idx), idx), normalqforward)];
                const int& normalqbackward =
                    TTraits::Stencil::QMap
                        .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                  data.getNeighbor(k, Stencil::Opposites[idx]))
                                  .NormalDirection)
                        ->second;
                double csolidbackward =
                    param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), normalqbackward)];

                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 *
                               (-(csolidforward) + 5 * param[data.getNeighbors()[k * Stencil::Q + idx]] - 3 * param[k] -
                                csolidbackward);

            } else {
                const int& normalq =
                    TTraits::Stencil::QMap
                        .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                  data.getNeighbor(data.getNeighbor(k, idx), idx))
                                  .NormalDirection)
                        ->second;
                double csolid = param[data.getNeighbor(data.getNeighbor(data.getNeighbor(k, idx), idx), normalq)];

                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 *
                               (-(csolid) + 5 * param[data.getNeighbors()[k * Stencil::Q + idx]] - 3 * param[k] -
                                param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]]);
            }

        } else {
            const int& normalq =
                TTraits::Stencil::QMap
                    .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                              data.getNeighbor(k, Stencil::Opposites[idx]))
                              .NormalDirection)
                    ->second;
            double csolid = param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), normalq)];

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 *
                           (-param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx]] +
                            5 * param[data.getNeighbors()[k * Stencil::Q + idx]] - 3 * param[k] - (csolid));
        }
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
}
