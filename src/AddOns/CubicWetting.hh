#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "AddOnBase.hh"
#include <math.h>

class CubicWetting : public AddOnBase {
    public:

        CubicWetting() = default;

        CubicWetting(const CubicWetting& other) : mPrefactor(other.mPrefactor) {}

        template<class TTraits>
        inline void compute(int k);

        template<class TTraits>
        inline void communicate();

        inline void setTheta(double theta);
        inline void setThetaDegrees(double theta);

    private:

        double mPrefactor = 0;
        

};

template<class TTraits>
inline void CubicWetting::compute(const int k) {

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using data = Data_Base<Lattice, Stencil>;

    const std::vector<int>& mv_Neighbors = data::getInstance().getNeighbors();

    if (Geometry<Lattice>::isSolid(k)) {

        double phiAvg = 0;
        int count = 0;

        for (int idx = 0; idx < Stencil::Q; idx++) {

            if (!Geometry<Lattice>::isSolid(mv_Neighbors[k * Stencil::Q+idx])) {

                phiAvg += OrderParameter<>::get<Lattice>(mv_Neighbors[k * Stencil::Q+idx]);
                count++;

            }

        }

        phiAvg /= count;
        OrderParameter<>::get<Lattice>(k) = phiAvg - mPrefactor * (pow(phiAvg,2) - 1.0);

    }
}


inline void CubicWetting::setTheta(const double theta){
    mPrefactor = cos(theta) / sqrt(2.0);
}


inline void CubicWetting::setThetaDegrees(const double theta){
    setTheta(theta / 180.0 * M_PI);
}


template<class TTraits>
inline void CubicWetting::communicate(){

    using Lattice = typename TTraits::Lattice;
    Lattice::communicate(OrderParameter<>::getInstance<Lattice>());
    
}
