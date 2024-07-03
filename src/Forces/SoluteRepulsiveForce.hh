#pragma once
#include <iostream>
#include <stdexcept>

#include "../Forcing.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "ForceBase.hh"

/**
 * \file SoluteRepulsiveForce.hh
 * \brief Contains the force class for an applied force based on the order parameter gradient.
 */

/**
 * \brief SoluteRepulsiveForce can be used to apply a force to a passive scalar dissolved in one phase
 * of a two-phase system. The force is proportional to the gradient of the order parameter and is necessary
 * to prevent the solute from diffusing into the other phase.
 */
template <class TMethod = Guo>
class SoluteRepulsiveForce : public ForceBase<TMethod> {
   public:
    //! Return force at lattice point k in direction xyz
    template <class TTraits>
    inline double computeXYZ(int xyz, int k);

    //! Return force at lattice point k along lattice direction idx
    template <class TTraits>
    inline double computeQ(int idx, int k);

   private:
    // Strength of the repulsive force
    int strength = 1.0;
};

template <class TMethod>
template <class TTraits>
inline double SoluteRepulsiveForce<TMethod>::computeXYZ(int xyz, int k) {
    using Lattice = typename TTraits::Lattice;

    double concentration = Solute<>::get<Lattice>(k);

    double gradOP_x = GradientOrderParameter<>::get<Lattice, Lattice::NDIM>(k, 0);
    double gradOP_y = GradientOrderParameter<>::get<Lattice, Lattice::NDIM>(k, 1);
    double gradOP_z = GradientOrderParameter<>::get<Lattice, Lattice::NDIM>(k, 2);

    if (xyz == 0) return gradOP_x * concentration * strength;
    if (xyz == 1) return gradOP_y * concentration * strength;
    if (xyz == 2) return gradOP_z * concentration * strength;
    throw std::invalid_argument("Invalid index for force component");
}

template <class TMethod>
template <class TTraits>
inline double SoluteRepulsiveForce<TMethod>::computeQ(int idx, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    double forceCi = computeXYZ(0, k) * Stencil::Ci_xyz(0)[idx];
    if constexpr (Lattice::NDIM > 1) forceCi += computeXYZ(1, k) * Stencil::Ci_xyz(1)[idx];
    if constexpr (Lattice::NDIM > 2) forceCi += computeXYZ(2, k) * Stencil::Ci_xyz(2)[idx];
    return forceCi;
}
