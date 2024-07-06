#pragma once
#pragma once
#include <math.h>

#include <iostream>
#include <utility>

#include "../Lattice.hh"
#include "../Parameters.hh"
#include "../Service.hh"
#include "ForceBase.hh"

// ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
// unfinished (should be able to specify magnitude and direction).

template <class TMethod = AllenCahnSourceMethod, int TComponentID = 0>
class AllenCahnSource : public ForceBase<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(int xyz, int k) const;  // Return force at TTraits::Lattice point k in direction xyz

    template <class TTraits>
    inline double computeQ(int xyz, int k) const;

    double mD;

    double mTau = 0.5;

    double mMobility = 0.00333;

    const double magnitudecutoff = 1e-14;

    inline void setD(double d) { mD = d; }
    inline void setTau(double tau) { mTau = (2 / tau - 1) / 2.0; }
    inline void setMobility(double mobility) { mMobility = mobility; }
    inline void setDTauAndMobility(double d, double tau, double mobility) {
        setD(d);
        setTau(tau);
        setMobility(mobility);
    }

    template <class TTraits>
    inline double computeBeta(int xyz, int k) const;
};

template <class TMethod, int TComponentID>
template <class TTraits>
inline double AllenCahnSource<TMethod, TComponentID>::computeXYZ(int xyz, int k) const {
    double gradx =
        GradientOrderParameter<TComponentID>::template get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0);
    double grady =
        GradientOrderParameter<TComponentID>::template get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 1);

    double magnitudegrad2 = gradx * gradx + grady * grady;
    if constexpr (TTraits::Lattice::NDIM == 3)
        magnitudegrad2 +=
            GradientOrderParameter<TComponentID>::template get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k,
                                                                                                                  2) *
            GradientOrderParameter<TComponentID>::template get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 2);
    double normal =
        GradientOrderParameter<TComponentID>::template get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) /
        sqrt(magnitudegrad2);

    if (sqrt(magnitudegrad2) > magnitudecutoff) {
        return (mMobility / mTau) *
               (4 * OrderParameter<TComponentID>::template get<typename TTraits::Lattice>(k) *
                    (1. - OrderParameter<TComponentID>::template get<typename TTraits::Lattice>(k)) * normal / mD -
                computeBeta<TTraits>(xyz, k));

    } else
        return 0;
}

template <class TMethod, int TComponentID>
template <class TTraits>
inline double AllenCahnSource<TMethod, TComponentID>::computeQ(int idx, int k) const {
    return 0;
}

template <class TMethod, int TComponentID>
template <class TTraits>
inline double AllenCahnSource<TMethod, TComponentID>::computeBeta(int xyz, int k) const {
    using Lattice = typename TTraits::Lattice;

    double sum = 0;
    double orderparametersum = 0;
    double gradientsum[TTraits::Lattice::NDIM] = {};
    for (int component = 0; component < TTraits::NumberOfComponents - 1; component++) {
        const double *gradOP =
            &getInstance<GradientOrderParameter, TTraits::NumberOfComponents - 1, Lattice, Lattice::NDIM>(
                component)[k * Lattice::NDIM];
        gradientsum[0] += gradOP[0];
        gradientsum[1] += gradOP[1];
        double magnitudegrad2 = pow(gradOP[0], 2) + pow(gradOP[1], 2);
        if constexpr (TTraits::Lattice::NDIM == 3) {
            magnitudegrad2 += pow(gradOP[2], 2);
            gradientsum[2] += gradOP[2];
        }
        double normal = gradOP[xyz] / sqrt(magnitudegrad2);
        if (sqrt(magnitudegrad2) > magnitudecutoff) {
            double orderParameter = getInstance<OrderParameter, TTraits::NumberOfComponents - 1, Lattice>(component)[k];
            orderparametersum += orderParameter;
            sum += 4 * orderParameter * (1 - orderParameter) * normal / mD;
        }
    }
    double magnitudegrad2 = gradientsum[0] * gradientsum[0] + gradientsum[1] * gradientsum[1];
    if constexpr (TTraits::Lattice::NDIM == 3) magnitudegrad2 += gradientsum[2] * gradientsum[2];
    double normal = (-gradientsum[xyz]) / (sqrt(magnitudegrad2));
    if (sqrt(magnitudegrad2) > magnitudecutoff) sum += 4 * (1 - orderparametersum) * (orderparametersum)*normal / mD;
    return OrderParameter<TComponentID>::template get<typename TTraits::Lattice>(k) * sum;
}
