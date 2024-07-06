#pragma once
#include <iostream>

#include "../Parameters.hh"
#include "BoundaryBase.hh"

class HaoCheng : public BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);

    template <class TTraits>
    inline void communicate(){};

    template <class TTraits, class TDistributionType>
    inline void communicate(TDistributionType& mDistribution);

    inline void setPhi(const double phi_in) { phi = phi_in; }

   private:
    // Default value for the order parameter at the boundary
    double phi = -1.0;
};

template <class TTraits, class TDistributionType>
inline void HaoCheng::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {
        double weightSum = 0.0;
        double gstar = phi;

        for (int q = 0; q < TTraits::Stencil::Q; q++) {
            if (!this->apply<Lattice>(distribution.streamIndex(k, q))) {
                weightSum += TTraits::Stencil::Weights[q];
            } else {
                gstar -= distribution.getDistributionPointer(k)[distribution.getOpposite(q)];
            }
        }

        if (!this->apply<Lattice>(distribution.streamIndex(k, idx))) {
            distribution.getDistributionPointer(k)[distribution.getOpposite(idx)] =
                gstar * TTraits::Stencil::Weights[distribution.getOpposite(idx)] / weightSum;
        }
    }
}

template <class TTraits, class TDistributionType>
inline void HaoCheng::communicate(TDistributionType& distribution) {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
}
