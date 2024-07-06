#pragma once
#include <iostream>

#include "../Parameters.hh"
#include "BoundaryBase.hh"

class VelocityInflow : public BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);

    template <class TTraits>
    inline void communicate(){};

    template <class TTraits, class TDistributionType>
    inline void communicate(TDistributionType& mDistribution);

    inline void setWallVelocity(const std::vector<double>& momentum) { mWallMomentum = momentum; };

   private:
    std::vector<double> mWallMomentum;
};

template <class TTraits, class TDistributionType>
inline void VelocityInflow::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {
        if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;

        double cidotmomentum = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            cidotmomentum += TTraits::Stencil::Ci_xyz(xyz)[idx] * mWallMomentum[xyz];
        }
        distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] =
            distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx), distribution.getOpposite(idx)) -
            2 * TTraits::Stencil::Weights[idx] * cidotmomentum / TTraits::Stencil::Cs2;
    }
}

template <class TTraits, class TDistributionType>
inline void VelocityInflow::communicate(TDistributionType& distribution) {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
}
