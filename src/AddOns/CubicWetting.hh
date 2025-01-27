#pragma once
#include <math.h>

#include <array>
#include <functional>
#include <map>

#include "../Geometry.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "AddOnBase.hh"
#include "Service.hh"

class CubicWetting : public AddOnBase {
   public:
    inline void setTheta(double theta);
    inline void setTheta(double (*theta)(int, int, int));

    inline void setThetaDegrees(double theta);
    inline void setThetaDegrees(double (*theta)(int, int, int));

    inline void setAlpha(double alpha);

    /**
     * \brief Set the thickness of the neutral wet layer at the inlet and outlet to avoid
     *       wetting phase towards them.
     * \param dir Direction of the neutral wet layer, typically main flow direction.
     *            Use x, y, or z to specify the direction.
     * \param thickness Thickness of the neutral wet layer in lattice units.
     *                  This value should be a non-negative integer.
     * \throws std::invalid_argument if the direction is not x, y, or z.
     */
    inline void setNeutralWetLayerThickness(const int dir, const int thickness) {
        if (dir != x && dir != y && dir != z) {
            throw std::invalid_argument("Invalid direction. Use x, y, or z.");
        }
        neutralWetLayerThickness = thickness;
        neutralWetLayerDirection = dir;
    }

    template <class TTraits>
    inline void compute(int k);

    template <class TTraits>
    inline void communicate();

    // Set the flag to true when simple or coupled mass loss models are used; othewise, unwanted phase change from
    // soli-fluid boundary will occur, resulting in the order parameter going beyond (-1, 1).
    bool useSinglePhaseCheck = false;

   private:
    double mAlpha = 2;
    double mPrefactor = 0;
    int neutralWetLayerThickness = 0;
    int neutralWetLayerDirection = 0;
    std::map<int, double> mPrefactorMap;
    std::function<double(std::array<int, 3>, double)> mPrefactorFn;
    enum { x = 0, y = 1, z = 2 };
};

template <class TTraits>
inline void CubicWetting::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    using data = Data_Base<Lattice, Stencil>;

    if (!this->apply<Lattice>(k)) return;

    // Calculate prefactor if using non-constant contact angle
    double prefactor;
    if (mPrefactorFn) {
        if (mPrefactorMap.find(k) == mPrefactorMap.end()) {
            auto xyz = computeXYZ<Lattice>(k);
            mPrefactorMap[k] = mPrefactorFn(xyz, mAlpha);
        }
        prefactor = mPrefactorMap[k];
    } else {
        prefactor = mPrefactor;
    }

    bool neighborsSinglePhase = useSinglePhaseCheck;

    // Get average order parameter from the neighbours
    double phiAvg = 0;
    int count = 0;
    for (int idx = 1; idx < Stencil::Q; idx++) {
        int neighbor = data::getInstance().getNeighbors()[k * Stencil::Q + idx];

        if (!Geometry<Lattice>::isBoundary(neighbor)) {
            phiAvg += OrderParameter<>::get<Lattice>(neighbor);
            count++;

            if (abs(OrderParameter<>::get<Lattice>(neighbor)) < 0.5) neighborsSinglePhase = false;
        }
    }
    phiAvg /= count;

    // Neutral wetting situations
    if (neighborsSinglePhase) {
        OrderParameter<>::get<Lattice>(k) = phiAvg;
        return;
    } else if (neutralWetLayerThickness > 0) {
        int coord, maxCoord;
        auto [coord_x, coord_y, coord_z] = computeXYZ<Lattice>(k);
        switch (neutralWetLayerDirection) {
            case x:
                coord = coord_x;
                maxCoord = Lattice::LX;
                break;
            case y:
                coord = coord_y;
                maxCoord = Lattice::LY;
                break;
            case z:
                coord = coord_z;
                maxCoord = Lattice::LZ;
                break;
            default:
                throw std::invalid_argument("Invalid neutralWetLayerDirection");
        }

        if (coord < neutralWetLayerThickness || coord >= maxCoord - neutralWetLayerThickness) {
            OrderParameter<>::get<Lattice>(k) = phiAvg;
            return;
        }
    }

    // Set the order parameter on the solid node
    OrderParameter<>::get<Lattice>(k) = phiAvg - prefactor * (pow(phiAvg, 2) - 1.0);
}

inline void CubicWetting::setTheta(double theta) { mPrefactor = cos(theta) / (sqrt(2.0) * mAlpha); }

inline void CubicWetting::setTheta(double (*theta)(int, int, int)) {
    mPrefactorFn = [theta](std::array<int, 3> xyz, double alpha) {
        double thetaK = theta(xyz[0], xyz[1], xyz[2]);
        return cos(thetaK) / (sqrt(2.0) * alpha);
    };
}

inline void CubicWetting::setThetaDegrees(double theta) { setTheta(theta / 180.0 * M_PI); }

inline void CubicWetting::setThetaDegrees(double (*theta)(int, int, int)) {
    mPrefactorFn = [theta](std::array<int, 3> xyz, double alpha) {
        double thetaK = theta(xyz[0], xyz[1], xyz[2]) / 180.0 * M_PI;
        return cos(thetaK) / (sqrt(2.0) * alpha);
    };
}

inline void CubicWetting::setAlpha(double alpha) {
    mPrefactor *= mAlpha / alpha;
    mAlpha = alpha;
}

template <class TTraits>
inline void CubicWetting::communicate() {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicate(OrderParameter<>::getInstance<Lattice>());
}
