#pragma once
#include <iostream>
#include <stdexcept>

#include "../Forcing.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "ForceBase.hh"

/**
 * \file ExternalForce.hh
 * \brief Contains the force class for a constant applied body force in a given direction.
 */

/**
 * \brief BodyForce can be used to apply an external force to the whole fluid or a single component
 */
template <class TMethod = Guo>
class BodyForce : public ForceBase<TMethod> {
   public:
    //! Set the force vector and, optionally, the fluid component to which it applies
    inline void setForce(std::vector<double> force, int component = -1);

    inline void setMagnitudeX(double magnitude);  //!< Set the x-component of the force
    inline void setMagnitudeY(double magnitude);  //!< Set the y-component of the force
    inline void setMagnitudeZ(double magnitude);  //!< Set the z-component of the force

    inline void activateGravityY();  //!< (DEPRECATED) Activates gravity in the y-direction.
    inline void activateGravityZ();  //!< (DEPRECATED) Activates gravity in the z-direction.

    //! Return force at lattice point k in direction xyz
    template <class TTraits>
    inline double computeXYZ(int xyz, int k);

    //! Return force at lattice point k along lattice direction idx
    template <class TTraits>
    inline double computeQ(int idx, int k);

    //! Calculate any possible source/correction term for velocity
    template <class TTraits>
    inline double computeVelocitySource(int xyz, int k);

   private:
    int mComponent = -1;
    double mMagnitudeX = 0;
    double mMagnitudeY = 0;
    double mMagnitudeZ = 0;
};

template <class TMethod>
inline void BodyForce<TMethod>::setForce(std::vector<double> force, int component) {
    mMagnitudeX = force[0];
    if (force.size() > 1) mMagnitudeY = force[1];
    if (force.size() > 2) mMagnitudeZ = force[2];
    mComponent = component;
}

template <class TMethod>
inline void BodyForce<TMethod>::setMagnitudeX(double magnitude) {
    mMagnitudeX = magnitude;
}

template <class TMethod>
inline void BodyForce<TMethod>::setMagnitudeY(double magnitude) {
    mMagnitudeY = magnitude;
}

template <class TMethod>
inline void BodyForce<TMethod>::setMagnitudeZ(double magnitude) {
    mMagnitudeZ = magnitude;
}

template <class TMethod>
inline void BodyForce<TMethod>::activateGravityY() {
    mComponent = 0;
    print("DEPRECATION WARNING: Use setForce on component 0 (or 1) for a gravitational force.");
}

template <class TMethod>
inline void BodyForce<TMethod>::activateGravityZ() {
    mComponent = 0;
    print("DEPRECATION WARNING: Use setForce on component 0 (or 1) for a gravitational force.");
}

template <class TMethod>
template <class TTraits>
inline double BodyForce<TMethod>::computeXYZ(int xyz, int k) {
    using Lattice = typename TTraits::Lattice;

    // Density of fluid being forced
    // TODO: Implement ternary / N-component
    double density = 1;
    if (mComponent == 0) {
        density = 0.5 + 0.5 * OrderParameter<>::get<Lattice>(k);
    } else if (mComponent == 1) {
        density = 0.5 - 0.5 * OrderParameter<>::get<Lattice>(k);
    }

    if (xyz == 0) return mMagnitudeX * density;
    if (xyz == 1) return mMagnitudeY * density;
    if (xyz == 2) return mMagnitudeZ * density;
    throw std::invalid_argument("Invalid index for force component");
}

template <class TMethod>
template <class TTraits>
inline double BodyForce<TMethod>::computeQ(int idx, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    double forceCi = computeXYZ(0, k) * Stencil::Ci_xyz(0)[idx];
    if constexpr (Lattice::NDIM > 1) forceCi += computeXYZ(1, k) * Stencil::Ci_xyz(1)[idx];
    if constexpr (Lattice::NDIM > 2) forceCi += computeXYZ(2, k) * Stencil::Ci_xyz(2)[idx];
    return forceCi;
}

template <class TMethod>
template <class TTraits>
inline double BodyForce<TMethod>::computeVelocitySource(int xyz, int k) {
    return +computeXYZ<TTraits>(xyz, k) * TTraits::Lattice::DT / (2.0);
}
