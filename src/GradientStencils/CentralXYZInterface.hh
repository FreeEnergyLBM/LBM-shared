#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZInterface : GradientBase<Gradient, Cartesian> {
    template <class TTraits, class TParameter>
    inline double compute(int direction, int k);

    inline void setInterfaceId(int id) { mInterfaceID = id; }

    inline void setInterfaceDistance(double (*distance)(int k, int idx)) { evalInterfaceDistance = distance; }

    inline void setInterfaceVal(double value) { mInterfaceVal = value; }

    double mInterfaceVal = 0;

    int mInterfaceID = 5;

    static double defaultDistance(int k, int idx) { return 0.5; }

    double (*evalInterfaceDistance)(int k, int idx) = &defaultDistance;
};

template <class TTraits, class TParameter>
inline double CentralXYZInterface::compute(int direction, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;

    if (Geometry<Lattice>::getBoundaryType(k) == 0) {
        for (int idx = 1; idx < Stencil::Q; idx++) {
            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) != mInterfaceID) &&
                (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx])) != mInterfaceID)) {
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] *
                               (TParameter::template get<Lattice>(data.getNeighbor(k, idx)));

            } else if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == mInterfaceID &&
                       (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx])) !=
                        mInterfaceID)) {
                double interfacedistance = evalInterfaceDistance(k, idx);
                gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (mInterfaceVal) /
                               (1 + interfacedistance);

            } else if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) != mInterfaceID &&
                       Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx])) ==
                           mInterfaceID) {
                double interfacedistance = evalInterfaceDistance(k, Stencil::Opposites[idx]);
                gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] *
                               (TParameter::template get<Lattice>(data.getNeighbor(k, idx))) / (1 + interfacedistance);
            }
        }
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
}
