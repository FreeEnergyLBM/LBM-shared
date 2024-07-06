#pragma once
#include <iostream>
#include <unordered_map>

#include "../Parameters.hh"
#include "BoundaryBase.hh"

class VelocityInflowVariable : public BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);

    template <class TTraits>
    inline void communicate(){};

    template <class TTraits, class TDistributionType>
    inline void communicate(TDistributionType& mDistribution);

    inline void setWallVelocity(const std::unordered_map<int, std::vector<double>>& momentum) {
        mWallMomentum = momentum;
    };

   private:
    std::unordered_map<int, std::vector<double>> mWallMomentum;
};

template <class TTraits, class TDistributionType>
inline void VelocityInflowVariable::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k) || !mWallMomentum.count(k)) return;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {
        if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;

        double cidotmomentum = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            cidotmomentum += TTraits::Stencil::Ci_xyz(xyz)[idx] * mWallMomentum[k][xyz];
        }

        distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] =
            distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx), distribution.getOpposite(idx)) -
            2 * TTraits::Stencil::Weights[idx] * cidotmomentum / TTraits::Stencil::Cs2;
    }
}

template <class TTraits, class TDistributionType>
inline void VelocityInflowVariable::communicate(TDistributionType& distribution) {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
}
