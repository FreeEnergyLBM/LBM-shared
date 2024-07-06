#pragma once
#include <iostream>

#include "../Collide.hh"
#include "../Parameters.hh"
#include "BoundaryBase.hh"

template <class TParam>
class Refill : public BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);

    template <class TTraits>
    inline void communicatePrecompute();

    template <class TTraits>
    inline void communicate(){};

    template <class TTraits, class TDistributionType>
    inline void communicate(TDistributionType& mDistribution);
};

template <class TParam>
template <class TTraits, class TDistributionType>
inline void Refill<TParam>::compute(TDistributionType& distribution, int k) {
    static_assert(TDistributionType::SaveEquilibrium,
                  "Cannot use this boundary type unless a data type that saves equilibrium information is used.");

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    using DataType = Data_Base<Lattice, Stencil>;
    DataType& data = DataType::getInstance();

    double avgparam = 0;
    int sum = 0;

    for (int idx = 0; idx < Stencil::Q; idx++) {
        int neighbor = data.getNeighbors()[k * Stencil::Q + idx];
        if (Geometry<Lattice>::getBoundaryType(neighbor) == 0) {
            avgparam += TParam::template get<Lattice>(neighbor);
            sum += 1;
        }
    }

    if (sum == 0) {
        for (int idx = 0; idx < Stencil::Q; idx++) {
            int neighbor = data.getNeighbors()[k * Stencil::Q + idx];
            if (Geometry<Lattice>::getBoundaryType(neighbor) == 0 ||
                Geometry<Lattice>::getBoundaryType(neighbor) == 6 ||
                Geometry<Lattice>::getBoundaryType(neighbor) == 5) {
                avgparam += TParam::template get<Lattice>(neighbor);
                sum += 1;
            }
        }
    }
    avgparam *= 1.0 / (double)sum;

    for (int idx = 0; idx < Stencil::Q; idx++) {
        double updatedist = Stencil::Weights[idx] * avgparam; /* *
                                 (1.0 + CollisionBase<Lattice, Stencil>::computeVelocityFactorFirstOrder(
                                            Velocity<>::getAddress<Lattice>(k, 0), idx)) +
                             (distribution.getDistributionPointer(k)[distribution.getOpposite(idx)] -
                              distribution.getEquilibriumPointer(k)[distribution.getOpposite(idx)]);*/

        distribution.getDistributionPointer(k)[idx] = updatedist;
        distribution.getDistributionOldPointer(k)[idx] = updatedist;
    }
}

template <class TParam>
template <class TTraits>
inline void Refill<TParam>::communicatePrecompute() {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicate(TParam::template getInstance<Lattice>());
}

template <class TParam>
template <class TTraits, class TDistributionType>
inline void Refill<TParam>::communicate(TDistributionType& distribution) {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
    Lattice::communicateDistributionAllEquilibrium(distribution);
}
