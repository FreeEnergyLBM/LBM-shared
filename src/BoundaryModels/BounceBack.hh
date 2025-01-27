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

// This class has the same functionality as BounceBack class, but applies the bounce-back BC by setting NodeID to fluid nodes.
class BounceBackWithFluidNodeID : public BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);
    inline void setBoundaryID(int id) { boundaryID = id; }

   private:
    int boundaryID = 1;
};

template <class TTraits, class TDistributionType>
inline void BounceBackWithFluidNodeID::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;

    for (int idx = 1; idx < TTraits::Stencil::Q; idx++) {
        if (Geometry<Lattice>::getBoundaryType(distribution.streamIndex(k, distribution.getOpposite(idx))) ==
            boundaryID) {
            distribution.getDistributionPointer(k)[idx] =
                distribution.getPostCollisionDistribution(k, distribution.getOpposite(idx));
        }
    }
}