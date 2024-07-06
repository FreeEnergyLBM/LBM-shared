#pragma once
#include <math.h>

#include <iostream>

#include "../Parameters.hh"
#include "AddOnBase.hh"

template <class TParameter,int TNdim=1>
class Communicate : public AddOnBase {
   public:

    template <class TTraits>
    inline void communicate();
};

template <class TParameter, int TNdim>
template <class TTraits>
inline void Communicate<TParameter,TNdim>::communicate() {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicate(TParameter::template getInstance<Lattice,TNdim>());
}