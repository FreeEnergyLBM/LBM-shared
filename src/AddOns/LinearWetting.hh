#pragma once
#include "../Parameters.hh"
#include "AddOnBase.hh"
#include<iostream>
#include<math.h>


class LinearWetting : public AddOnBase {
    public:

        LinearWetting() = default;

        LinearWetting(const LinearWetting& other) : m_Theta(other.m_Theta), m_Omega(other.m_Omega), m_Prefactor(other.m_Prefactor) {}

        template<typename traits>
        inline void compute(const int k);

        template<typename traits>
        inline void communicate();

    private:

        inline double calcOmega(const double theta);

        double m_Theta = M_PI / 2.0;
        double m_Omega=0;
        double m_Prefactor = 0.0;

    public:
        inline void setTheta(const double theta);

        inline void setThetaDegrees(const double theta);

        inline void setPrefactor(const double prefactor);

        inline void setPrefactor(const double A, const double kappa);

};

template<typename traits>
inline void LinearWetting::compute(const int k) {

    Geometry<typename traits::Lattice> m_Geometry;

    using data = Data_Base<typename traits::Lattice, typename traits::Stencil>;

    if (m_Geometry.isSolid(k)) {
        
        double wettingsum = 0;
        int count = 0;

        for (int idx = 0; idx < traits::Stencil::Q; idx++) {
                
            if (!m_Geometry.isSolid(data::getInstance().getNeighbors()[k * traits::Stencil::Q+idx])) {

                wettingsum += m_Prefactor * m_Omega + OrderParameter<>::get<typename traits::Lattice>(data::getInstance().getNeighbors()[k * traits::Stencil::Q+idx]);
                count++;

            }

        }

        OrderParameter<>::get<typename traits::Lattice>(k) = wettingsum/((double)count);

    }
    
}

inline void LinearWetting::setTheta(const double theta){

    m_Theta = theta;
    m_Omega = calcOmega(m_Theta);

}

inline void LinearWetting::setThetaDegrees(const double theta){

    m_Theta = M_PI * theta / 180.0;
    m_Omega = calcOmega(m_Theta);

}

inline double LinearWetting::calcOmega(const double theta){

    double alpha = acos(sin(theta) * sin(theta));
    return 2 * (((M_PI / 2.0 - theta) >= 0) - ((M_PI / 2.0 - theta) < 0)) * sqrt(cos(alpha / 3.0) * (1.0 - cos(alpha / 3.0)));

}

inline void LinearWetting::setPrefactor(const double prefactor){

    m_Prefactor=prefactor;

}

inline void LinearWetting::setPrefactor(const double A, const double kappa){

    m_Prefactor=sqrt(A/(2.0*kappa));

}

template<typename traits>
inline void LinearWetting::communicate(){

    using Lattice = typename traits::Lattice;
    Lattice::communicate(OrderParameter<>::getInstance<Lattice>());
    
}
