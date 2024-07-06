#pragma once
#include <math.h>

#include <iostream>
#include <utility>

#include "../Geometry.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "../Service.hh"
#include "AddOnBase.hh"

template <class TParam, class TOrderParameter>
class CalculateFromOrderParameter : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

    inline void setValues(std::vector<double> values) { mValues = values; }

    std::vector<double> mValues = {};
};

template <class TTraits>
inline void CalculateFromOrderParameter::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    const int N = TTraits::NumberOfComponents;

    if (Geometry<Lattice>::isBulkSolid(k)) return;

    double sum = mValues.back();

    for (int component = 0; component < N - 1; component++) {
        double chemPot = getInstance<ChemicalPotential, N - 1, Lattice>(component)[k];
        double gradOP = getGradientInstance<TGradientType, OrderParameter, N - 1, Lattice, TDirections>(
            component)[k * TDirections + idx];
        sum += getInstance<TOrderParameter, N, Lattice>(component)[k] * (mValues[component] - mValues.back());
    }

    TParam::get<Lattice>(k) = 1.0 / sum;
}
