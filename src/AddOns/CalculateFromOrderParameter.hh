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

    if (Geometry<Lattice>::isBulkSolid(k)) return;

    double sum = mValues.back();

    for (int component = 0; component < TTraits::NumberOfComponents - 1; component++) {
        // std::cout<<ChemicalPotential<>::template get<typename TTraits::Lattice>(k)<<std::endl;
        double chemPot =
            ChemicalPotential<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice>(k, component);
        double gradOP = TGradientType<OrderParameter<TTraits::NumberOfComponents - 1>,
                                      (TTraits::NumberOfComponents - 1)>::template get<typename TTraits::Lattice,
                                                                                       TDirections>(k, component, idx);
        sum += TOrderParameter::get<Lattice>(k, component) * (mValues[component] - mValues.back());
        //
    }

    TParam::get<Lattice>(k) = 1.0 / sum;
}