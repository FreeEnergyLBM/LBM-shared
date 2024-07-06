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

    // TODO by Mehrdad: This should be written in a more general way to allow for inclusion of y and z directions if
    // needed!
    /**
     * \brief Set the thickness of the neutral wet layer at the inlet and outlet to avoid
     *       wetting phase towards them.
     * \param thickness Thickness of the neutral wet layer.
     */
    inline void setNeutralWetLayerThickness(const int thickness) { neutralWetLayerThickness = thickness; }

    template <class TTraits>
    inline void compute(int k);

    template <class TTraits>
    inline void communicate();

    bool useSinglePhaseCheck = false;

   private:
    double mAlpha = 2;
    double mPrefactor = 0;
    int neutralWetLayerThickness = 0;
    std::map<int, double> mPrefactorMap;
    std::function<double(std::array<int, 3>, double)> mPrefactorFn;
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
        int x = computeXGlobal<Lattice>(k);
        if (x < neutralWetLayerThickness || x >= Lattice::LX - neutralWetLayerThickness) {
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
