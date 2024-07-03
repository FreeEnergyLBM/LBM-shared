#pragma once
#pragma once
#include <math.h>

#include <iostream>
#include <utility>

#include "../Lattice.hh"
#include "../Parameters.hh"
#include "../Service.hh"
#include "ForceBase.hh"

// Implementation of the AD sink/source term according to sect 8.3.5.1, Kruger book

template <class TMethod>
class ADSinkSource : public ForceBase<TMethod> {
   public:
    template <class TTraits>
    inline double compute(int k) const;

    template <class TTraits>
    inline double computeDensitySource(int k) const;

    inline void setSinkSourceTerm(double rate) { mSourceRate = rate; }

   private:
    double mSourceRate = 0.0;
};

template <class TMethod>
template <class TTraits>
inline double ADSinkSource<TMethod>::compute(int k) const {
    return mSourceRate;
}

template <class TMethod>
template <class TTraits>
inline double ADSinkSource<TMethod>::computeDensitySource(int k) const {
    return +0.5 * TTraits::Lattice::DT * mSourceRate;
}
