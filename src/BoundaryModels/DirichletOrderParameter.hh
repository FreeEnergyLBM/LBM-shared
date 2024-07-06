#pragma once
#include <iostream>

#include "../Parameters.hh"
#include "BoundaryBase.hh"

class DirichletOrderParameter : public BoundaryBase {
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
inline void DirichletOrderParameter::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {
        if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;

        distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] =
            -distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),
                                                       distribution.getOpposite(idx)) +
            2 * TTraits::Stencil::Weights[idx] * ChemicalPotential<>::get<Lattice>(distribution.streamIndex(k, idx)) /
                TTraits::Stencil::Cs2 * mInterfaceVal;
    }
}

template <class TTraits, class TDistributionType>
inline void DirichletOrderParameter::communicate(TDistributionType& distribution) {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
}
