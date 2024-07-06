#pragma once
#pragma once
#include <iostream>
#include <utility>

#include "../Forcing.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "../Service.hh"
#include "ForceBase.hh"

// ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
// unfinished (should be able to specify magnitude and direction).

template <class TMethod = Guo, template <class> class TGradientType = Gradient>
class ChemicalForceBinary : public ForceBase<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(int xyz, int k);  // Return force at traits::Lattice point k in direction xyz

    template <class TTraits>
    inline double computeQ(int xyz, int k);

    template <class TTraits>
    inline double computeVelocitySource(int xyz, int k);  // Calculate any possible source/correction term for velocity

    template <class TTraits, int TDirections>
    inline double computeChemicalForce(int idx, int k);
};

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceBinary<TMethod, TGradientType>::computeXYZ(int xyz, int k) {
    if constexpr (has_type<Cartesian, typename TMethod::mt_Stencils>::type::value) {
        return computeChemicalForce<TTraits, TTraits::Lattice::NDIM>(xyz, k);
    }
    return 0;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceBinary<TMethod, TGradientType>::computeQ(int idx, int k) {
    if constexpr (has_type<AllDirections, typename TMethod::mt_Stencils>::type::value) {
        return computeChemicalForce<TTraits, TTraits::Stencil::Q>(idx, k);
    }
    return 0;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits, int TDirections>
inline double ChemicalForceBinary<TMethod, TGradientType>::computeChemicalForce(int idx, int k) {
    using Lattice = typename TTraits::Lattice;
    constexpr int N = TTraits::NumberOfComponents;

    double sum = 0;
    for (int component = 0; component < N - 1; component++) {
        double chemPot = getInstance<ChemicalPotential, N - 1, Lattice>(component)[k];
        double gradOP = getGradientInstance<TGradientType, OrderParameter, N - 1, Lattice, TDirections>(
            component)[k * TDirections + idx];
        sum += chemPot * gradOP;
    }

    return sum;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceBinary<TMethod, TGradientType>::computeVelocitySource(
    const int xyz, const int k) {  // Need to correct velocity

    return +computeXYZ<TTraits>(xyz, k) * TTraits::Lattice::DT / (2.0);
}

template <class TMethod = Guo, template <class> class TGradientType = Gradient>
class ChemicalForceBinaryMu : public ForceBase<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(int xyz, int k);  // Return force at traits::Lattice point k in direction xyz

    template <class TTraits>
    inline double computeQ(int xyz, int k);

    template <class TTraits>
    inline double computeVelocitySource(int xyz, int k);  // Calculate any possible source/correction term for velocity

    template <class TTraits, int TDirections>
    inline double computeChemicalForce(int idx, int k);

    template <class TTraits>
    inline void communicate();
};

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceBinaryMu<TMethod, TGradientType>::computeXYZ(int xyz, int k) {
    if constexpr (has_type<Cartesian, typename TMethod::mt_Stencils>::type::value) {
        return computeChemicalForce<TTraits, TTraits::Lattice::NDIM>(xyz, k);
    }
    return 0;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceBinaryMu<TMethod, TGradientType>::computeQ(int idx, int k) {
    if constexpr (has_type<AllDirections, typename TMethod::mt_Stencils>::type::value) {
        return computeChemicalForce<TTraits, TTraits::Stencil::Q>(idx, k);
    }
    return 0;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits, int TDirections>
inline double ChemicalForceBinaryMu<TMethod, TGradientType>::computeChemicalForce(int idx, int k) {
    using Lattice = typename TTraits::Lattice;
    constexpr int N = TTraits::NumberOfComponents;

    double sum = 0;
    for (int component = 0; component < N - 1; component++) {
        double orderParam = getInstance<OrderParameter, N - 1, Lattice>(component)[k];
        double gradChemPot = getGradientInstance<TGradientType, ChemicalPotential, N, Lattice, TDirections>(
            component)[k * TDirections + idx];
        sum += orderParam * gradChemPot;
    }

    return -sum;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceBinaryMu<TMethod, TGradientType>::computeVelocitySource(
    const int xyz, const int k) {  // Need to correct velocity

    return +computeXYZ<TTraits>(xyz, k) * TTraits::Lattice::DT / (2.0);
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline void ChemicalForceBinaryMu<TMethod, TGradientType>::communicate() {  // Not necessary

    using Lattice = typename TTraits::Lattice;
    Lattice::communicate(ChemicalPotential<>::getInstance<Lattice>());
}

template <class TMethod = Guo, template <class> class TGradientType = Gradient>
class ChemicalForce : public ForceBase<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(int xyz, int k);  // Return force at traits::Lattice point k in direction xyz

    template <class TTraits>
    inline double computeQ(int xyz, int k);

    template <class TTraits>
    inline double computeVelocitySource(int xyz, int k);  // Calculate any possible source/correction term for velocity

    template <class TTraits, int TDirections>
    inline double computeChemicalForce(int idx, int k);
};

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForce<TMethod, TGradientType>::computeXYZ(int xyz, int k) {
    if constexpr (has_type<Cartesian, typename TMethod::mt_Stencils>::type::value) {
        return computeChemicalForce<TTraits, TTraits::Lattice::NDIM>(xyz, k);
    }
    return 0;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForce<TMethod, TGradientType>::computeQ(int idx, int k) {
    if constexpr (has_type<AllDirections, typename TMethod::mt_Stencils>::type::value) {
        return computeChemicalForce<TTraits, TTraits::Stencil::Q>(idx, k);
    }
    return 0;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits, int TDirections>
inline double ChemicalForce<TMethod, TGradientType>::computeChemicalForce(int idx, int k) {
    using Lattice = typename TTraits::Lattice;
    constexpr int N = TTraits::NumberOfComponents;

    double sum = 0;
    double gradopsum = 0;

    for (int component = 0; component < TTraits::NumberOfComponents - 1; component++) {
        // TMP_INSTANCE: Previously the total number of instances was N-(N==2), why? is this an issue?
        double chemPot = getInstance<ChemicalPotential, N, Lattice>(component)[k];
        double gradOP = getGradientInstance<TGradientType, OrderParameter, N - 1, Lattice, TDirections>(
            component)[k * TDirections + idx];
        sum += chemPot * gradOP;
        gradopsum += gradOP;
    }

    if constexpr ((TTraits::NumberOfComponents > 2))
        sum += ChemicalPotential<N - 1>::template get<Lattice>(k) * (-gradopsum);

    return sum;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForce<TMethod, TGradientType>::computeVelocitySource(const int xyz,
                                                                           const int k) {  // Need to correct velocity

    return +computeXYZ<TTraits>(xyz, k) * TTraits::Lattice::DT / (2.0);
}

template <class TMethod = Guo, template <class> class TGradientType = Gradient>
class ChemicalForceMu : public ForceBase<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(int xyz, int k);  // Return force at traits::Lattice point k in direction xyz

    template <class TTraits>
    inline double computeQ(int xyz, int k);

    template <class TTraits>
    inline double computeVelocitySource(int xyz, int k);  // Calculate any possible source/correction term for velocity

    template <class TTraits, int TDirections>
    inline double computeChemicalForce(int idx, int k);
};

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceMu<TMethod, TGradientType>::computeXYZ(int xyz, int k) {
    if constexpr (has_type<Cartesian, typename TMethod::mt_Stencils>::type::value) {
        return computeChemicalForce<TTraits, TTraits::Lattice::NDIM>(xyz, k);
    }
    return 0;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceMu<TMethod, TGradientType>::computeQ(int idx, int k) {
    if constexpr (has_type<AllDirections, typename TMethod::mt_Stencils>::type::value) {
        return computeChemicalForce<TTraits, TTraits::Stencil::Q>(idx, k);
    }
    return 0;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits, int TDirections>
inline double ChemicalForceMu<TMethod, TGradientType>::computeChemicalForce(int idx, int k) {
    using Lattice = typename TTraits::Lattice;
    constexpr int N = TTraits::NumberOfComponents;

    double sum = 0;
    double opsum = 0;

    for (int component = 0; component < N - 1; component++) {
        double OP = getInstance<OrderParameter, N, Lattice>(component)[k];
        double gradchemPot = getGradientInstance<TGradientType, ChemicalPotential, N, Lattice, TDirections>(
            component)[k * TDirections + idx];
        sum += OP * gradchemPot;
        opsum += OP;
    }

    if constexpr ((TTraits::NumberOfComponents > 2))
        sum += TGradientType<ChemicalPotential<N - 1>>::template get<Lattice, TDirections>(k, idx) * (1.0 - opsum);

    return sum;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceMu<TMethod, TGradientType>::computeVelocitySource(const int xyz,
                                                                             const int k) {  // Need to correct velocity

    return +computeXYZ<TTraits>(xyz, k) * TTraits::Lattice::DT / (2.0);
}

template <class TMethod = Guo, template <class> class TGradientType = Gradient>
class ChemicalForceRho : public ForceBase<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(int xyz, int k);  // Return force at traits::Lattice point k in direction xyz

    template <class TTraits>
    inline double computeQ(int xyz, int k);

    template <class TTraits>
    inline double computeVelocitySource(int xyz, int k);  // Calculate any possible source/correction term for velocity

    template <class TTraits, int TDirections>
    inline double computeChemicalForce(int idx, int k);
};

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceRho<TMethod, TGradientType>::computeXYZ(int xyz, int k) {
    if constexpr (has_type<Cartesian, typename TMethod::mt_Stencils>::type::value) {
        return computeChemicalForce<TTraits, TTraits::Lattice::NDIM>(xyz, k);
    }
    return 0;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceRho<TMethod, TGradientType>::computeQ(int idx, int k) {
    if constexpr (has_type<AllDirections, typename TMethod::mt_Stencils>::type::value) {
        return computeChemicalForce<TTraits, TTraits::Stencil::Q>(idx, k);
    }
    return 0;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits, int TDirections>
inline double ChemicalForceRho<TMethod, TGradientType>::computeChemicalForce(int idx, int k) {
    using Lattice = typename TTraits::Lattice;
    constexpr int N = TTraits::NumberOfComponents;

    double sum = 0;
    for (int component = 0; component < N - 1; component++) {
        double chemPot = getInstance<ChemicalPotential, N - 1, Lattice>(component)[k];
        double gradDensity =
            getGradientInstance<TGradientType, Density, N - 1, Lattice>(component)[k * Lattice::NDIM + idx];
        sum += chemPot * gradDensity;
    }

    return sum;
}

template <class TMethod, template <class> class TGradientType>
template <class TTraits>
inline double ChemicalForceRho<TMethod, TGradientType>::computeVelocitySource(
    const int xyz, const int k) {  // Need to correct velocity

    return +computeXYZ<TTraits>(xyz, k) * TTraits::Lattice::DT / (2.0);
}
