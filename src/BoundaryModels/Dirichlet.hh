#pragma once
#include <iostream>

#include "../Parameters.hh"
#include "BoundaryBase.hh"

class Dirichlet : public BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);

    template <class TTraits>
    inline void communicate(){};

    template <class TTraits, class TDistributionType>
    inline void communicate(TDistributionType& mDistribution);

    inline void setInterfaceVal(double val) { mInterfaceVal = val; };

   private:
    double mInterfaceVal;
};

template <class TTraits, class TDistributionType>
inline void Dirichlet::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;

    for (int idx = 1; idx < TTraits::Stencil::Q; idx++) {
        if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;

        distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] =
            -distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),
                                                       distribution.getOpposite(idx)) +
            2 * TTraits::Stencil::Weights[idx] * mInterfaceVal;
    }
}

template <class TTraits, class TDistributionType>
inline void Dirichlet::communicate(TDistributionType& distribution) {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
}
