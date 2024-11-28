#pragma once
#include <iostream>
#include <stdexcept>

#include "../Forcing.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "ForceBase.hh"

/**
 * \file CapillaryForcePorous.hh
 * @brief CapillaryForcePorous class is to compute the capillary force in porous media.
 * The capillary force is responsible for the movement of the weting phase, and is computed
 * as the gradient of the capillary pressure.
 */

template <class TMethod = GuoPorous, template <class> class TGradientType = Gradient,
          class TParameterScalar = WettingRelativePermeability<>, class TParameterVector = Velocity<>>
class CapillaryForcePorous : public ForceBase<TMethod> {
   public:
    //! Return force at lattice point k in direction xyz
    template <class TTraits>
    inline double computeXYZ(int xyz, int k);

    //! Return force at lattice point k along lattice direction idx
    template <class TTraits>
    inline double computeQ(int idx, int k);

    //! Calculate any possible source/correction term for velocity
    template <class TTraits>
    inline double computeVelocitySource(int xyz, int k);

    template <class TTraits, int TDirections>
    inline double computeCapillaryForce(int idx, int k);

   private:
    int mComponent = -1;
    double darcyPrefactor = 0;
    double forchheimerPrefactor = 0;

    enum { x = 0, y = 1, z = 2 };  // Indices corresponding to x, y, z directions
};

template <class TMethod, template <class> class TGradientType, class TParameterScalar, class TParameterVector>
template <class TTraits>
inline double CapillaryForcePorous<TMethod, TGradientType, TParameterScalar, TParameterVector>::computeXYZ(int xyz,
                                                                                                           int k) {
    using Lattice = typename TTraits::Lattice;

    std::vector<double>& relativePermeability =
        TParameterScalar::template get<Lattice>();  // Reference to vector of relative permeabilities

    std::vector<double>& velocity =
        TParameterVector::template get<Lattice, Lattice::NDIM>();  // Reference to vector of velocities

    double saturation = Saturation<>::get<typename TTraits::Lattice>(k);
    double porosity = Porosity<>::get<typename TTraits::Lattice>(k);
    double permeability = Permeability<>::get<typename TTraits::Lattice>(k);

    auto lbModel = static_cast<ModelBase<typename TTraits::Lattice, TTraits>*>(this->mModel);
    double kinematicViscosity = TTraits::Stencil::Cs2 * (lbModel->getTau() - 0.5 * TTraits::Lattice::DT);

    std::vector localVelocity{0.0, 0.0, 0.0};
    localVelocity[x] = velocity[k * TTraits::Stencil::D + x];
    if constexpr (TTraits::Lattice::NDIM > 1) localVelocity[y] = velocity[k * TTraits::Stencil::D + y];
    if constexpr (TTraits::Lattice::NDIM > 2) localVelocity[z] = velocity[k * TTraits::Stencil::D + z];

    double localVelocityMagnitude = sqrt(localVelocity[x] * localVelocity[x] + localVelocity[y] * localVelocity[y] +
                                         localVelocity[z] * localVelocity[z]);

    // capillary force term
    double gradPc = 0;
    if constexpr (has_type<Cartesian, typename TMethod::mt_Stencils>::type::value)
        gradPc = computeCapillaryForce<TTraits, TTraits::Lattice::NDIM>(xyz, k);
    double capillaryForce = 1 / saturation * gradPc;

    // TODO: SEPARATE THE DARCYS AND FORCHHEIMERS FORCES.
    if (localVelocityMagnitude != 0) {
        // Darcy Prefactor
        darcyPrefactor = -porosity * kinematicViscosity / (permeability * relativePermeability[k]);

        // Forchheimer Prefactor
        double F_eps = 1.75 / sqrt(150. * porosity * porosity * porosity);
        forchheimerPrefactor =
            -porosity * F_eps * localVelocityMagnitude / sqrt(permeability * relativePermeability[k]);
    }

    return darcyPrefactor * localVelocity[xyz] + forchheimerPrefactor * localVelocity[xyz] + porosity * capillaryForce;
    throw std::invalid_argument("Invalid index for force component");
}

template <class TMethod, template <class> class TGradientType, class TParameterScalar, class TParameterVector>
template <class TTraits>
inline double CapillaryForcePorous<TMethod, TGradientType, TParameterScalar, TParameterVector>::computeQ(int idx,
                                                                                                         int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    double forceCi = computeXYZ(0, k) * Stencil::Ci_xyz(0)[idx];
    if constexpr (Lattice::NDIM > 1) forceCi += computeXYZ(1, k) * Stencil::Ci_xyz(1)[idx];
    if constexpr (Lattice::NDIM > 2) forceCi += computeXYZ(2, k) * Stencil::Ci_xyz(2)[idx];
    return forceCi;
}

template <class TMethod, template <class> class TGradientType, class TParameterScalar, class TParameterVector>
template <class TTraits, int TDirections>
inline double CapillaryForcePorous<TMethod, TGradientType, TParameterScalar, TParameterVector>::computeCapillaryForce(
    int idx, int k) {
    using Lattice = typename TTraits::Lattice;
    constexpr int N = TTraits::NumberOfComponents;

    double sum = 0;
    for (int component = 0; component < N - 1; component++) {
        double gradPc = getGradientInstance<TGradientType, CapillaryPressure, N - 1, Lattice, TDirections>(
            component)[k * TDirections + idx];
        sum += gradPc;
    }

    return sum;
}

template <class TMethod, template <class> class TGradientType, class TParameterScalar, class TParameterVector>
template <class TTraits>
inline double CapillaryForcePorous<TMethod, TGradientType, TParameterScalar, TParameterVector>::computeVelocitySource(
    int xyz, int k) {
    return +computeXYZ<TTraits>(xyz, k) * TTraits::Lattice::DT / (2.0);
}
