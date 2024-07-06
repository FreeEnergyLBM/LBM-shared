#pragma once
#include <iostream>

#include "BoundaryBase.hh"

class BounceBack : public BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);
};

template <class TTraits, class TDistributionType>
inline void BounceBack::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;

    for (int idx = 1; idx < TTraits::Stencil::Q; idx++) {
        if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;

        distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] =
            distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx), distribution.getOpposite(idx));
    }
}
