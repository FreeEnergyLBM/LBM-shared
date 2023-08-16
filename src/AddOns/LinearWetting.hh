#pragma once
#include "../Parameters.hh"
#include "../Geometry.hh"
#include "AddOnBase.hh"
#include<iostream>
#include<math.h>

class LinearWetting : public AddOnBase {
    public:

        LinearWetting() = default;

        LinearWetting(const LinearWetting& other) : mTheta(other.mTheta), mOmega(other.mOmega), mPrefactor(other.mPrefactor) {}

        template<class TTraits>
        inline void compute(int k);

        template<class TTraits>
        inline void communicate();

    private:

        inline double calcOmega(double theta);

        double mTheta = M_PI / 2.0;
        double mOmega=0;
        double mPrefactor = 0.0;

    public:
        inline void setTheta(double theta);

        inline void setThetaDegrees(double theta);

        inline void setPrefactor(double prefactor);

        inline void setPrefactor(double A, double kappa);

};

template<class TTraits>
inline void LinearWetting::compute(int k) {

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using data = Data_Base<Lattice, Stencil>;

    if (Geometry<Lattice>::getBoundaryType(k)==1) {
        
        double wettingsum = 0;
        int count = 0;

        for (int idx = 0; idx < Stencil::Q; idx++) {
                
            if (Geometry<Lattice>::getBoundaryType(data::getInstance().getNeighbors()[k * Stencil::Q + idx])!=1) {

                wettingsum += mPrefactor * mOmega + OrderParameter<>::get<Lattice>(data::getInstance().getNeighbors()[k * Stencil::Q + idx]);
                count++;

            }

        }

        OrderParameter<>::get<Lattice>(k) = wettingsum / ((double)count);

    }
    
}

inline void LinearWetting::setTheta(double theta){

    mTheta = theta;
    mOmega = calcOmega(mTheta);

}

inline void LinearWetting::setThetaDegrees(double theta){

    mTheta = M_PI * theta / 180.0;
    mOmega = calcOmega(mTheta);

}

inline double LinearWetting::calcOmega(double theta){

    double alpha = acos(sin(theta) * sin(theta));
    return 2 * (((M_PI / 2.0 - theta) >= 0) - ((M_PI / 2.0 - theta) < 0)) * sqrt(cos(alpha / 3.0) * (1.0 - cos(alpha / 3.0)));

}

inline void LinearWetting::setPrefactor(double prefactor){

    mPrefactor = prefactor;

}

inline void LinearWetting::setPrefactor(double A, double kappa){

    mPrefactor = sqrt(A / (2.0 * kappa));

}

template<class TTraits>
inline void LinearWetting::communicate(){

    using Lattice = typename TTraits::Lattice;
    Lattice::communicate(OrderParameter<>::getInstance<Lattice>());
    
}
