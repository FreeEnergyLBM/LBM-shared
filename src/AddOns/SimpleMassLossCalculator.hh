#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "AddOnBase.hh"
#include <math.h>

/**
 * @brief A class for evaporation with a constant rate per unit area.
 */
class SimpleMassLossCalculator : public AddOnBase {
    public:
        /// Compute the volumetric evaporation rate
        template<class TTraits>
        inline void compute(int k);

        /// Set the constant evaporation rate per unit area
        inline void setEvaporationRate(double rate);

        inline void setFluidPhase(int phase) {mFluidPhase=phase;}

    private:
        double mEvaporationRate = 1e-4;
        int mFluidPhase=0;
};


template<class TTraits>
inline void SimpleMassLossCalculator::compute(int k){
    if (Geometry<typename TTraits::Lattice>::isBulkSolid(k)) return;

    // Get the local gradient of the order parameter
    double gradOP = 0;
    for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
        gradOP += pow(GradientOrderParameter<TTraits::NumberOfComponents-1>::template get<typename TTraits::Lattice,TTraits::Lattice::NDIM>(k, mFluidPhase, xyz), 2);
    }
    gradOP = sqrt(gradOP);

    // Calculate the volumetric mass loss
    MassSink<>::get<typename TTraits::Lattice>(k) = mEvaporationRate * gradOP;
}


inline void SimpleMassLossCalculator::setEvaporationRate(double rate){
    mEvaporationRate = rate;
}