#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "AddOnBase.hh"
#include <math.h>

class CubicWetting : public AddOnBase {
    public:

        CubicWetting() = default;

        CubicWetting(const CubicWetting& other) : m_Prefactor(other.m_Prefactor) {}

        template<typename traits>
        inline void compute(const int k);

        template<typename traits>
        inline void communicate();

        inline void setTheta(const double theta);
        inline void setThetaDegrees(const double theta);

    private:

        double m_Prefactor = 0;
        

};

template<typename traits>
inline void CubicWetting::compute(const int k) {
    using data = Data_Base<typename traits::Lattice, typename traits::Stencil>;
    const std::vector<int>& mv_Neighbors = data::getInstance().getNeighbors();
    Geometry<typename traits::Lattice> m_Geometry;
    if (m_Geometry.isSolid(k)) {
        double phiAvg = 0;
        int count = 0;
        for (int idx = 0; idx < traits::Stencil::Q; idx++) {
            if (!m_Geometry.isSolid(mv_Neighbors[k * traits::Stencil::Q+idx])) {
                phiAvg += OrderParameter<>::get<typename traits::Lattice>(mv_Neighbors[k * traits::Stencil::Q+idx]);
                count++;
            }
        }
        phiAvg /= count;
        OrderParameter<>::get<typename traits::Lattice>(k) = phiAvg - m_Prefactor * (pow(phiAvg,2) - 1.0);
    }
}


inline void CubicWetting::setTheta(const double theta){
    m_Prefactor = cos(theta) / sqrt(2.0);
}


inline void CubicWetting::setThetaDegrees(const double theta){
    setTheta(theta / 180.0 * M_PI);
}


template<typename traits>
inline void CubicWetting::communicate(){

    using Lattice = typename traits::Lattice;
    Lattice::communicate(OrderParameter<>::getInstance<Lattice>());
    
}
