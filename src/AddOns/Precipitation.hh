#pragma once
#include <iostream>

#include "AddOnBase.hh"

template <class TParameter>
class Precipitation : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

    inline void activatePrecipitation() { mPrecipitationActive = true; }
    inline void setLowPrecipitationThreshold(double value) { mPrecipitationThresholdLow = value; }
    inline void setHighPrecipitationThreshold(double value) { mPrecipitationThresholdHigh = value; }

   private:
    double mPrecipitationThresholdLow = 0.8;
    double mPrecipitationThresholdHigh = 1.0;
    bool mPrecipitationActive = true;
};

template <class TParameter>
template <class TTraits>
inline void Precipitation<TParameter>::compute(int k) {
    using Lattice = typename TTraits::Lattice;

    if (Geometry<Lattice>::getBoundaryType(k) == 0) {
        if (mPrecipitationActive) {
            if (TParameter::template get<Lattice>(k) > mPrecipitationThresholdLow &&
                TParameter::template get<Lattice>(k) < mPrecipitationThresholdHigh) {
                Geometry<Lattice>::modifyBoundaryLabel(k, 1);
                // std::cout << "Precipitation at node " << k << "\n";
            }
        }
    }
}