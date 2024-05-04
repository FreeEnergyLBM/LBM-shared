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

    inline void setLiquidID(int id) { mLiquidID = id; }
    inline void setGasID(int id) { mGasID = id; }

   private:
    double mEvaporationRate = 1e-4;
    int mLiquidID = 0;
    int mGasID = 0;
};

template <class TTraits>
inline void SimpleMassLossCalculator::compute(int k) {
    if (Geometry<typename TTraits::Lattice>::isBulkSolid(k)) return;

    // Get the local gradient of the order parameter
    double gradOP = 0;
    if (mGasID < TTraits::NumberOfComponents - 1) {
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            gradOP +=
                fabs(GradientOrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice,
                                                                                           TTraits::Lattice::NDIM>(
                         k, mLiquidID, xyz) *
                     GradientOrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice,
                                                                                           TTraits::Lattice::NDIM>(
                         k, mGasID, xyz));
        }
    } else {
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            gradOP += fabs(
                GradientOrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice,
                                                                                      TTraits::Lattice::NDIM>(
                    k, mLiquidID, xyz) *
                (GradientOrderParameter<TTraits::NumberOfComponents -
                                        1>::template get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0, xyz) +
                 GradientOrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice,
                                                                                       TTraits::Lattice::NDIM>(k, 1,
                                                                                                               xyz)));
        }
    }

    gradOP = sqrt(gradOP);

    // Calculate the volumetric mass loss
    MassSink<>::get<typename TTraits::Lattice>(k) = mEvaporationRate * gradOP;
}

inline void SimpleMassLossCalculator::setEvaporationRate(double rate) { mEvaporationRate = rate; }