#pragma once
#include <iostream>

#include "../Parameters.hh"
#include "BoundaryBase.hh"

class InterpolatedDirichlet : public BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);

    template <class TTraits>
    inline void communicate(){};

    template <class TTraits, class TDistributionType>
    inline void communicate(TDistributionType& mDistribution);

    inline void setInterfaceVal(double val) { mInterfaceVal = val; };

    inline void setInterfaceDistanceFunction(double (*func)(int idx, int k)) { evalDistanceFunction = func; }

   private:
    static double defaultDistanceFunction(int idx, int k) { return true; }

    double (*evalDistanceFunction)(int idx, int k) = &defaultDistanceFunction;

    double mInterfaceVal;
};

template <class TTraits, class TDistributionType>
inline void InterpolatedDirichlet::compute(TDistributionType& distribution, int k) {  // Modify to include wall velocity

    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;

    for (int idx = 1; idx < TTraits::Stencil::Q; idx++) {
        if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;

        double dist = evalDistanceFunction(distribution.streamIndex(k, idx), distribution.getOpposite(idx));

        if (dist <= 0.5) {
            if (Geometry<Lattice>::isBoundary(distribution.streamIndex(distribution.streamIndex(k, idx), idx))) {
                distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] =
                    -distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),
                                                               distribution.getOpposite(idx)) +
                    2 * TTraits::Stencil::Weights[idx] * mInterfaceVal;
            } else {
                distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] =
                    -2 * (dist)*distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),
                                                                          distribution.getOpposite(idx)) +
                    (2 * dist - 1) * distribution.getPostCollisionDistribution(
                                         distribution.streamIndex(distribution.streamIndex(k, idx), idx),
                                         distribution.getOpposite(idx)) +
                    2 * TTraits::Stencil::Weights[idx] * mInterfaceVal;
            }
        } else
            distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] =
                -1.0 / (2.0 * dist) *
                    distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),
                                                              distribution.getOpposite(idx)) +
                (1 - 1.0 / (2.0 * dist)) *
                    distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx), idx) +
                2 * (1.0 / (2.0 * dist)) * TTraits::Stencil::Weights[idx] * mInterfaceVal;
    }
}

template <class TTraits, class TDistributionType>
inline void InterpolatedDirichlet::communicate(TDistributionType& distribution) {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
    Lattice::communicate(GradientOrderParameter<>::getInstance<Lattice, Lattice::NDIM>());
    Lattice::communicate(BoundaryLabels<Lattice::NDIM>::template getInstance<Lattice>());
}
