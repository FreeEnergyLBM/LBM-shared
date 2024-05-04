#pragma once
#include <math.h>

#include <iostream>

#include "../Parameters.hh"
#include "AddOnBase.hh"

template <class TParameter, class TParameterOld>
class SetParameterOld : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

    template <class TTraits>
    inline void communicate();
};

template <class TParameter, class TParameterOld>
template <class TTraits>
inline void SetParameterOld<TParameter, TParameterOld>::compute(int k) {
    using Lattice = typename TTraits::Lattice;

    TParameterOld::template get<Lattice>(k) = TParameter::template get<Lattice>(k);
}

template <class TParameter, class TParameterOld>
template <class TTraits>
inline void SetParameterOld<TParameter, TParameterOld>::communicate() {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicate(TParameterOld::template getInstance<Lattice>());
}
