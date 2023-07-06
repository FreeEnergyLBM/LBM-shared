#pragma once
#include "../Parameters.hh"
#include "AddOnBase.hh"
#include<iostream>
#include<math.h>


class LinearWetting : public AddOnBase {
    public:

        LinearWetting() = default;

        LinearWetting(const LinearWetting& other) : m_Theta(other.m_Theta), m_Omega(other.m_Omega), m_Prefactor(other.m_Prefactor) {}

        template<typename T_traits>
        inline void compute(int k);

        template<typename T_traits>
        inline void communicate();

    private:

        inline double calcOmega(double theta);

        double m_Theta = M_PI / 2.0;
        double m_Omega=0;
        double m_Prefactor = 0.0;

    public:
        inline void setTheta(double theta);

        inline void setThetaDegrees(double theta);

        inline void setPrefactor(double prefactor);

        inline void setPrefactor(double A, double kappa);

};

template<typename T_traits>
inline void LinearWetting::compute(int k) {

    using Lattice = typename T_traits::Lattice;
    using Stencil = typename T_traits::Stencil;

    using data = Data_Base<Lattice, Stencil>;

    if (Geometry<Lattice>::isSolid(k)) {
        
        double wettingsum = 0;
        int count = 0;

        for (int idx = 0; idx < Stencil::Q; idx++) {
                
            if (!Geometry<Lattice>::isSolid(data::getInstance().getNeighbors()[k * Stencil::Q + idx])) {

                wettingsum += m_Prefactor * m_Omega + OrderParameter<>::get<Lattice>(data::getInstance().getNeighbors()[k * Stencil::Q + idx]);
                count++;

            }

        }

        OrderParameter<>::get<Lattice>(k) = wettingsum / ((double)count);

    }
    
}

inline void LinearWetting::setTheta(double theta){

    m_Theta = theta;
    m_Omega = calcOmega(m_Theta);

}

inline void LinearWetting::setThetaDegrees(double theta){

    m_Theta = M_PI * theta / 180.0;
    m_Omega = calcOmega(m_Theta);

}

inline double LinearWetting::calcOmega(double theta){

    double alpha = acos(sin(theta) * sin(theta));
    return 2 * (((M_PI / 2.0 - theta) >= 0) - ((M_PI / 2.0 - theta) < 0)) * sqrt(cos(alpha / 3.0) * (1.0 - cos(alpha / 3.0)));

}

inline void LinearWetting::setPrefactor(double prefactor){

    m_Prefactor = prefactor;

}

inline void LinearWetting::setPrefactor(double A, double kappa){

    m_Prefactor = sqrt(A / (2.0 * kappa));

}

template<typename T_traits>
inline void LinearWetting::communicate(){

    using Lattice = typename T_traits::Lattice;
    Lattice::communicate(OrderParameter<>::getInstance<Lattice>());
    
}
