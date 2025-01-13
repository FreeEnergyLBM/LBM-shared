#pragma once
#include <math.h>

#include "../Lattice.hh"
#include "../Parameters.hh"
#include "AddOnBase.hh"

/**
 * @brief This class calculates the mass exchange between FELBM lattice and HLBM lattice.
 */

/**
 * @brief A class for evaporation coupled with the vapour phase evolution.
 */
class MassExchangeCalculator : public AddOnBase {
   public:
    /// Compute the volumetric evaporation rate
    template <class TTraits>
    inline void compute(int k);

    inline void setSwIr(double Sw_ir) { this->Sw_ir = Sw_ir; }

   private:
    double Sw_ir = 0.4;
};

template <class TTraits>
inline void MassExchangeCalculator::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    if (Geometry<Lattice>::isBulkSolid(k)) return;

    using DataType = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;
    DataType& data = DataType::getInstance();

    // if condition to ensure that the processor is only applied to salt crystal
    if (Geometry<Lattice>::getBoundaryType(k) == 7) {
        // loop over Q directions to find the nodes with boundary type == 0
        // and the amount to be exchanged.

        // Liquid mass of nodes with boundary type == 7
        double currentmass = Saturation<>::get<typename TTraits::Lattice>(k);

        // Liquid exchangaable mass of nodes with boundary type == 0
        int numberofzeronodes = 0;
        double localexchangepotential = 0.0;
        for (int idx = 1; idx < TTraits::Stencil::Q; idx++) {
            if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == 0) {
                numberofzeronodes++;
                localexchangepotential +=
                    (1.0 + OrderParameter<>::get<typename TTraits::Lattice>(data.getNeighbor(k, idx))) / 2.0 -
                    currentmass;
            }
        }

        // Do this only if there is any node with boundary type == 0
        if (numberofzeronodes > 0) {
            {
                // For cases that localexchangepotential > 0:
                double massexchange = 0.0;
                if (localexchangepotential > 0)
                    // Choose the min of localexchangepotential and (1.0 - currentmass)
                    massexchange = std::min(localexchangepotential, 1.0 - currentmass);

                // For cases that localexchangepotential < 0:
                else
                    // Choose the max of localexchangepotential and S_wi - currentmass
                    massexchange = std::max(localexchangepotential, Sw_ir - currentmass);

                // Update the mass exchange for the node with boundary type == 7
                MassExchange<>::get<typename TTraits::Lattice>(k) = massexchange;

                // Update the mass exchange for the nodes with boundary type == 0
                for (int idx = 1; idx < TTraits::Stencil::Q; idx++) {
                    if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == 0) {
                        double weight = 0.0;
                        if (localexchangepotential != 0)
                            weight =
                                ((1.0 + OrderParameter<>::get<typename TTraits::Lattice>(data.getNeighbor(k, idx))) /
                                     2.0 -
                                 currentmass) /
                                localexchangepotential;
                        MassExchange<>::get<typename TTraits::Lattice>(data.getNeighbor(k, idx)) -=
                            massexchange * weight;
                    }
                }
            }
        }
    }
}

class ResetMassExchangeParameter : public AddOnBase {
   public:
    /// Compute the volumetric evaporation rate
    template <class TTraits>
    inline void compute(int k);
};

template <class TTraits>
inline void ResetMassExchangeParameter::compute(int k) {
    
    // The MassExchange parameter should be reset to zero at the end of 
    // each iteration to avoid accumulation of values.
    MassExchange<>::get<typename TTraits::Lattice>(k) = 0.0;
}
