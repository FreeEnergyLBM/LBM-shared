#pragma once
#include <math.h>

#include <iostream>
#include <utility>

#include "../Lattice.hh"
#include "../Parameters.hh"
#include "../Service.hh"
#include "ForceBase.hh"

// Implementation of the AD sink/source term according to sect 8.3.5.1, Kruger book coupled with evaporation

template <class TMethod>
class CoupledADSinkSource : public ForceBase<TMethod> {
   public:
    template <class TTraits>
    inline double compute(int k) const;

    template <class TTraits>
    inline double computeDensitySource(int k) const;
};

template <class TMethod>
template <class TTraits>
inline double CoupledADSinkSource<TMethod>::compute(int k) const {
    return MassSink<>::get<typename TTraits::Lattice>(k);
}

template <class TMethod>
template <class TTraits>
inline double CoupledADSinkSource<TMethod>::computeDensitySource(int k) const {
    return +0.5 * TTraits::Lattice::DT * MassSink<>::get<typename TTraits::Lattice>(k);
}
