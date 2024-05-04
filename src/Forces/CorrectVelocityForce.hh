#pragma once
#include <iostream>

#include "../Forcing.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "ForceBase.hh"

// ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
// unfinished (should be able to specify magnitude and direction).

template <class TMethod = SimpleForcing>
class CorrectVelocityForce : public ForceBase<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(int xyz, int k);  // Return force at lattice point k in direction xyz

   private:
    std::vector<double> mWallVelocity = {0, 0, 0};
    int mSolidPhase = 1;
};

template <class TMethod>
template <class TTraits>
inline double CorrectVelocityForce<TMethod>::computeXYZ(int xyz, int k) {
    using Lattice = typename TTraits::Lattice;

    if (xyz == 0)
        return OrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice>(k,
                                                                                                        mSolidPhase) *
               (mWallVelocity[0] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0)) *
               Density<>::get<typename TTraits::Lattice>(k) / TTraits::Lattice::DT;
    if constexpr (Lattice::NDIM == 2) {
        return OrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice>(k,
                                                                                                        mSolidPhase) *
               (mWallVelocity[1] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 1)) *
               Density<>::get<typename TTraits::Lattice>(k) / TTraits::Lattice::DT;
    }

    else if constexpr (Lattice::NDIM == 3) {
        if (xyz == 1)
            return OrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice>(
                       k, mSolidPhase) *
                   (mWallVelocity[1] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 1)) *
                   Density<>::get<typename TTraits::Lattice>(k) / TTraits::Lattice::DT;
        return OrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice>(k,
                                                                                                        mSolidPhase) *
               (mWallVelocity[2] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 2)) *
               Density<>::get<typename TTraits::Lattice>(k) / TTraits::Lattice::DT;
    }

    return OrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice>(k, mSolidPhase) *
           (mWallVelocity[0] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0)) *
           Density<>::get<typename TTraits::Lattice>(k) / TTraits::Lattice::DT;
}
