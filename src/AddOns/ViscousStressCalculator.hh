#pragma once
#include "../Geometry.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "../Service.hh"
#include "AddOnBase.hh"

class ViscousStressCalculator : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

};

template <class TTraits>
inline void ViscousStressCalculator::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    constexpr int numstresses = Lattice::NDIM * (Lattice::NDIM + 1) / 2;

    if (Geometry<Lattice>::isBulkSolid(k)) return;

    for (int dir = 0; dir < numstresses; dir++) {
        ViscousStress<>::get<Lattice, numstresses>(k, dir) = 0;
    }

    double tauloc = 1.0 / InverseTau<>::get<Lattice>(k);

    for (int idx = 0; idx < Stencil::Q; idx++) {

        for (int dir1 = 0; dir1 < Lattice::NDIM; dir1++) {
            for (int dir2 = 0; dir2 <= dir1; dir2++) {
                ViscousStress<>::get<Lattice, numstresses>(k, (dir1 + 1) * (dir2 + 1) - 1) =
                    (GradientVelocity<>::get<Lattice,Lattice::NDIM,Lattice::NDIM>(k,dir1,dir2)-GradientVelocity<>::get<Lattice,Lattice::NDIM,Lattice::NDIM>(k,dir2,dir1)); // Order will be xx, xy, yy, xz, yz, zz
            }
        }
    }

    double sum = 0;
    for (int dir1 = 0; dir1 < Lattice::NDIM; dir1++) {
        for (int dir2 = 0; dir2 <= dir1; dir2++) {
            ViscousStress<>::get<Lattice, numstresses>(k, (dir1 + 1) * (dir2 + 1) - 1) *= Density<>::get<Lattice>(k)*Stencil::Cs2*(tauloc-Lattice::DT*0.5);
            sum += pow(ViscousStress<>::get<Lattice, numstresses>(k, (dir1 + 1) * (dir2 + 1) - 1),2) *
                       (2 - 1 * (dir1 == dir2));
        }
    }

    ViscousDissipation<>::get<Lattice>(k) = (1.0 / (2.0 * Density<>::get<Lattice>(k) * (tauloc - 0.5) / 3.0)) * sum;
}