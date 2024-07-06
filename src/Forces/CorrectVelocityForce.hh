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
    constexpr int NCOMP = TTraits::NumberOfComponents;

    double solidOP = getInstance<OrderParameter, NCOMP, Lattice>(mSolidPhase)[k];
    if (xyz == 0)
        return solidOP * (mWallVelocity[0] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0)) *
               Density<>::get<typename TTraits::Lattice>(k) / TTraits::Lattice::DT;
    if constexpr (Lattice::NDIM == 2) {
        return solidOP * (mWallVelocity[1] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 1)) *
               Density<>::get<typename TTraits::Lattice>(k) / TTraits::Lattice::DT;
    }

    else if constexpr (Lattice::NDIM == 3) {
        if (xyz == 1)
            return solidOP *
                   (mWallVelocity[1] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 1)) *
                   Density<>::get<typename TTraits::Lattice>(k) / TTraits::Lattice::DT;
        return solidOP * (mWallVelocity[2] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 2)) *
               Density<>::get<typename TTraits::Lattice>(k) / TTraits::Lattice::DT;
    }

    return solidOP * (mWallVelocity[0] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0)) *
           Density<>::get<typename TTraits::Lattice>(k) / TTraits::Lattice::DT;
}
