#pragma once
#include <math.h>

#include "../Parameters.hh"
#include "AddOnBase.hh"

template <class TParameter,int TNDir=1>
class MirrorBoundary : public AddOnBase {
   public:
    MirrorBoundary() { this->setNodeID(4); }

    template <class TTraits>
    inline void compute(int k);
};

template <class TParameter,int TNDir>
template <class TTraits>
inline void MirrorBoundary<TParameter,TNDir>::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    const int NDIM = Lattice::NDIM;

    using DataType = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    if (!this->apply<Lattice>(k)) return;

    const std::vector<int>& neighbors = DataType::getInstance().getNeighbors();
    const std::array<int8_t, NDIM>& normal = BoundaryLabels<NDIM>::template get<Lattice>(k).NormalDirection;
    int idx = Stencil::QMap.find(normal)->second;
    if (idx > 0) {
        double magnormal = 0;

        for (int xyz = 0; xyz < Lattice::NDIM; xyz++) {
            magnormal += normal[xyz] * normal[xyz];
        }

        magnormal = sqrt(magnormal);
        for (int dir = 0; dir < TNDir; dir++) {
            if (Geometry<Lattice>::isCorner(k)) {
                for (int idx1 = 1; idx1 < Stencil::Q; idx1++) {
                    double cidotnormal = 0;
                    for (int xyz = 0; xyz < Lattice::NDIM; xyz++) {
                        cidotnormal += Stencil::Ci_xyz(xyz)[idx1] * normal[xyz];
                    }
                    if (cidotnormal > 0) {
                        TParameter::template get<Lattice,TNDir>(k,dir) =
                            TParameter::template get<Lattice,TNDir>(neighbors[k * Stencil::Q + idx1],dir);
                        TParameter::template get<Lattice,TNDir>(neighbors[k * Stencil::Q + Stencil::Opposites[idx1]],dir) = TParameter::template get<Lattice>(
                                neighbors[neighbors[k * Stencil::Q + idx1] * Stencil::Q + idx1],dir);
                    }
                }
            } else {
                TParameter::template get<Lattice,TNDir>(k,dir) = TParameter::template get<Lattice,TNDir>(neighbors[k * Stencil::Q + idx],dir);

                TParameter::template get<Lattice,TNDir>(data.getNeighbor(k, Stencil::Opposites[idx]),dir) =
                    TParameter::template get<Lattice,TNDir>(neighbors[neighbors[k * Stencil::Q + idx] * Stencil::Q + idx],dir);
            }
        }
    }
}
