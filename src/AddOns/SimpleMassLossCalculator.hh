#pragma once
#include <math.h>

#include "../Lattice.hh"
#include "../Parameters.hh"
#include "AddOnBase.hh"

/**
 * @brief A class for evaporation with a constant rate per unit area.
 */
class SimpleMassLossCalculator : public AddOnBase {
   public:
    /// Compute the volumetric evaporation rate
    template <class TTraits>
    inline void compute(int k);

    /// Set the constant evaporation rate per unit area
    inline void setEvaporationRate(double rate);

    // Set the saturation humidity
    inline void setSaturationHumidity(double humidity) { humidityMax = humidity; }

    inline void setLiquidID(int id) { mLiquidId = id; }
    inline void setGasID(int id) { mGasId = id; }

   private:
    double mEvaporationRate = 0.0;
    int mLiquidId = 0;
    int mGasId = 0;
    double humidityMax = 1.0;
};

template <class TTraits>
inline void SimpleMassLossCalculator::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    if (!this->apply<Lattice>(k)) return;

    constexpr int N = TTraits::NumberOfComponents;
    const std::vector<double>& gradLiquid = getInstance<GradientOrderParameter, N, Lattice, Lattice::NDIM>(mLiquidId);

    // Get the local gradient of the order parameter
    double gradOP = 0;
    if (mGasId < TTraits::NumberOfComponents - 1) {
        const std::vector<double>& gradGas = getInstance<GradientOrderParameter, N, Lattice, Lattice::NDIM>(mGasId);
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            gradOP += fabs(gradLiquid[k * Lattice::NDIM + xyz] * gradGas[k * Lattice::NDIM + xyz]);
        }
    } else {
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            gradOP += fabs(
                gradLiquid[k * Lattice::NDIM + xyz] *
                (GradientOrderParameter<0>::template get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) +
                 GradientOrderParameter<1>::template get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz)));
        }
    }

    gradOP = sqrt(gradOP);

    // Calculate the volumetric mass loss
    MassSink<>::get<typename TTraits::Lattice>(k) = mEvaporationRate * gradOP;
}

inline void SimpleMassLossCalculator::setEvaporationRate(double rate) { mEvaporationRate = rate; }

/**
 * @brief A class for evaporation coupled with the vapour phase evolution.
 */
class CoupledMassLossCalculator : public AddOnBase {
   public:
    /// Compute the volumetric evaporation rate
    template <class TTraits>
    inline void compute(int k);

    /// Set the constant evaporation rate per unit area
    inline void setEvaporationRate(double rate);

    // Set the saturation humidity
    inline void setSaturationHumidity(double humidity) { humidityMax = humidity; }

   private:
    double mEvaporationRate = 0.0;
    int mLiquidId = 0;
    int mGasId = 0;
    double humidityMax = 1.0;
};

template <class TTraits>
inline void CoupledMassLossCalculator::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    if (!this->apply<Lattice>(k)) return;

    constexpr int N = TTraits::NumberOfComponents;
    const std::vector<double>& gradLiquid = getInstance<GradientOrderParameter, N, Lattice, Lattice::NDIM>(mLiquidId);

    // Get the local gradient of the order parameter
    double gradOP = 0;
    if (mGasId < TTraits::NumberOfComponents - 1) {
        const std::vector<double>& gradGas = getInstance<GradientOrderParameter, N, Lattice, Lattice::NDIM>(mGasId);
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            gradOP += fabs(gradLiquid[k * Lattice::NDIM + xyz] * gradGas[k * Lattice::NDIM + xyz]);
        }
    } else {
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            gradOP += fabs(
                gradLiquid[k * Lattice::NDIM + xyz] *
                (GradientOrderParameter<0>::template get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) +
                 GradientOrderParameter<1>::template get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz)));
        }
    }

    gradOP = sqrt(gradOP);

    double localHumidity = Humidity<>::get<typename TTraits::Lattice>(k);

    // liquid mass per unit volume
    // NOTE: c1 is not the liquid content. It is exactly the amount of liquid in a lattice point.
    // So, it has the unit of LM/LL3.
    // Note that we assume \rho = 1.0 here.
    // double c1 = (1.0 + OrderParameter<>::get<typename TTraits::Lattice>(k)) / 2.0;

    // Calculate the volumetric mass loss
    MassSink<>::get<typename TTraits::Lattice>(k) =
        mEvaporationRate * 1.0 * (gradOP / 1.0) * (humidityMax - localHumidity);
}

inline void CoupledMassLossCalculator::setEvaporationRate(double rate) { mEvaporationRate = rate; }
