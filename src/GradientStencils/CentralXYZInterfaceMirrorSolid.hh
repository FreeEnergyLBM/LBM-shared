#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZInterfaceMirrorSolid : GradientBase<Gradient, Cartesian> {
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
inline double CentralXYZInterfaceMirrorSolid::compute(int direction, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;
    const static auto& param = TParameter::template get<Lattice>();

    if (Geometry<Lattice>::getBoundaryType(k) == 0 || Geometry<Lattice>::getBoundaryType(k) == 6) {
        for (int idx = 1; idx < Stencil::Q; idx++) {
            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == 1 ||
                 Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == 7)) {
                if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx])) == 0 ||
                    Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx])) == 4) {
                    const int& normalq =
                        TTraits::Stencil::QMap
                            .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                      data.getNeighbor(k, idx))
                                      .NormalDirection)
                            ->second;

                    double csolid = param[data.getNeighbor(data.getNeighbor(k, idx), normalq)];

                    gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * csolid;

                } else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx])) == 1))
                    ;
                else {
                    const int& normalq =
                        TTraits::Stencil::QMap
                            .find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(
                                      data.getNeighbor(k, idx))
                                      .NormalDirection)
                            ->second;

                    double csolid = param[data.getNeighbor(data.getNeighbor(k, idx), normalq)];

                    double interfacedistance = evalInterfaceDistance(k, Stencil::Opposites[idx]);
                    gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (csolid) /
                                   (1 + interfacedistance);
                }

            } else {
                if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == 0 ||
                     Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == 6 ||
                     Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == 4) &&
                    (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx])) !=
                     mInterfaceID)) {
                    gradientsum +=
                        Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (param[data.getNeighbor(k, idx)]);

                } else if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == mInterfaceID &&
                           (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx])) == 0 ||
                            Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx])) == 6)) {
                    double interfacedistance = evalInterfaceDistance(k, idx);
                    gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (mInterfaceVal) /
                                   (1 + interfacedistance);

                } else if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) != mInterfaceID &&
                           Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx])) ==
                               mInterfaceID) {
                    double interfacedistance = evalInterfaceDistance(k, Stencil::Opposites[idx]);
                    gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] *
                                   (param[data.getNeighbor(k, idx)]) / (1 + interfacedistance);

                } else if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == mInterfaceID &&
                           Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx])) == 1) {
                    double interfacedistance = evalInterfaceDistance(k, idx);
                    gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (mInterfaceVal) /
                                   (1 + interfacedistance);
                }
            }
        }
    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
}
