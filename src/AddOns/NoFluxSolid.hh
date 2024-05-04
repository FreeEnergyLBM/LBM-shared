#pragma once
#include <iostream>

#include "AddOnBase.hh"

template <class TParam>
class NoFluxSolid : public AddOnBase {
   public:
    NoFluxSolid() { this->setNodeID(5); }

    template <class TTraits>
    inline void compute(int k);
};

template <class TParam>
template <class TTraits>
inline void NoFluxSolid<TParam>::compute(int k) {
    if (Geometry<typename TTraits::Lattice>::isBulkSolid(k)) return;

    using Lattice = typename TTraits::Lattice;
    using DataType = Data_Base<Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    if (!this->apply<Lattice>(k)) return;

    double val = 0;
    double count = 0;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {
        int kNeighbor1 = data.getNeighbor(k, idx);
        int kNeighbor2 = data.getNeighbor(kNeighbor1, idx);
        if (this->apply<Lattice>(kNeighbor2)) {
            continue;
        } else if (this->apply<Lattice>(kNeighbor1)) {
            val += TParam::template get<Lattice>(kNeighbor2);
        } else {
            val += TParam::template get<Lattice>(kNeighbor1);
        }
        count += 1.0;
    }

    if (count > 0) TParam::template get<Lattice>(k) = val / count;
}
