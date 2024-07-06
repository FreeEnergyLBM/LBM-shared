#pragma once
#include <iostream>

#include "../Parameters.hh"
#include "BoundaryBase.hh"

class ExtrapolationOutflow : public BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);

    template <class TTraits>
    inline void communicate(){};

    template <class TTraits, class TDistributionType>
    inline void communicate(TDistributionType& mDistribution);
};

template <class TTraits, class TDistributionType>
inline void ExtrapolationOutflow::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    const int& normalq =
        Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;

        double cidotnormal = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            cidotnormal +=
                TTraits::Stencil::Ci_xyz(xyz)[idx] *
                BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection[xyz];
        }

        if (cidotnormal != 0) {
            distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] =
                (4.0 * distribution.getDistributionPointer(
                           distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[idx] -
                 distribution.getDistributionPointer(distribution.streamIndex(
                     distribution.streamIndex(distribution.streamIndex(k, normalq), normalq), normalq))[idx]) /
                3.0;

        } else {
            distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] =
                distribution.getEquilibriumPointer(distribution.streamIndex(k, normalq))[idx];
        }
    }
}

template <class TTraits, class TDistributionType>
inline void ExtrapolationOutflow::communicate(TDistributionType& distribution) {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
    Lattice::communicateDistributionAllEquilibrium(distribution);
}
