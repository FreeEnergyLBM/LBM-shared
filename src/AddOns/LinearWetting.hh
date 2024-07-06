#pragma once
#include <math.h>

#include <array>
#include <functional>
#include <map>

#include "../Geometry.hh"
#include "../Parameters.hh"
#include "AddOnBase.hh"
#include "Service.hh"

class LinearWetting : public AddOnBase {
   public:
    inline void setTheta(double theta);
    inline void setTheta(double (*theta)(int, int, int));

    inline void setThetaDegrees(double theta);
    inline void setThetaDegrees(double (*theta)(int, int, int));

    inline void setPrefactor(double prefactor);
    inline void setPrefactor(double A, double kappa);

    template <class TTraits>
    inline void compute(int k);

    template <class TTraits>
    inline void communicate();

   private:
    static inline double calcOmega(double theta);

    double mPrefactor = 0.0;
    double mOmega = 0;
    std::map<int, double> mOmegaMap;
    std::function<double(std::array<int, 3>)> mOmegaFn;
};

template <class TTraits>
inline void LinearWetting::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    using data = Data_Base<Lattice, Stencil>;

    if (!this->apply<Lattice>(k)) return;

    // Calculate omega if using non-constant contact angle
    double omega;
    if (mOmegaFn) {
        if (mOmegaMap.find(k) == mOmegaMap.end()) {
            auto xyz = computeXYZ<Lattice>(k);
            mOmegaMap[k] = mOmegaFn(xyz);
        }
        omega = mOmegaMap[k];
    } else {
        omega = mOmega;
    }

    // Get average order parameter from the neighbours
    double phiAvg = 0;
    int count = 0;
    for (int idx = 0; idx < Stencil::Q; idx++) {
        int neighbor = data::getInstance().getNeighbors()[k * Stencil::Q + idx];

        if (!Geometry<Lattice>::isBoundary(neighbor)) {
            phiAvg += OrderParameter<>::get<Lattice>(neighbor);
            count++;
        }
    }
    phiAvg /= count;

    // Set the order parameter on the solid node
    OrderParameter<>::get<Lattice>(k) = phiAvg + mPrefactor * omega;
}

inline void LinearWetting::setTheta(double theta) { mOmega = calcOmega(theta); }

inline void LinearWetting::setTheta(double (*theta)(int, int, int)) {
    mOmegaFn = [theta](std::array<int, 3> xyz) {
        double thetaK = theta(xyz[0], xyz[1], xyz[2]);
        return LinearWetting::calcOmega(thetaK);
    };
}

inline void LinearWetting::setThetaDegrees(double theta) { mOmega = calcOmega(theta / 180.0 * M_PI); }

inline void LinearWetting::setThetaDegrees(double (*theta)(int, int, int)) {
    mOmegaFn = [theta](std::array<int, 3> xyz) {
        double thetaK = theta(xyz[0], xyz[1], xyz[2]) / 180.0 * M_PI;
        return LinearWetting::calcOmega(thetaK);
    };
}

inline double LinearWetting::calcOmega(double theta) {
    double alpha = acos(sin(theta) * sin(theta));
    return 2 * (((M_PI / 2.0 - theta) >= 0) - ((M_PI / 2.0 - theta) < 0)) *
           sqrt(cos(alpha / 3.0) * (1.0 - cos(alpha / 3.0)));
}

inline void LinearWetting::setPrefactor(double prefactor) { mPrefactor = prefactor; }

inline void LinearWetting::setPrefactor(double A, double kappa) { mPrefactor = sqrt(A / (2.0 * kappa)); }

template <class TTraits>
inline void LinearWetting::communicate() {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicate(OrderParameter<>::getInstance<Lattice>());
}
