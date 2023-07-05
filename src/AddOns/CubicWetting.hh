#pragma once
#include "../Parameters.hh"
#include "AddOnBase.hh"
#include <math.h>



template<typename TLattice, typename TStencil>
class CubicWetting : public AddOnBase {
    public:

        inline void postprocess(const int k);

        inline void setTheta(const double theta);
        inline void setThetaDegrees(const double theta);

        inline void communicatePostProcess();

    private:

        OrderParameter<TLattice> m_OrderParameter;
        Geometry<TLattice> m_Geometry;
        Data_Base<TLattice, TStencil> m_Data;

        double m_Prefactor = 0;
        const std::vector<int> mv_Neighbors = m_Data.getNeighbors();

};


template<typename TLattice, typename TStencil>
inline void CubicWetting<TLattice,TStencil>::postprocess(const int k) {
    if (m_Geometry.isSolid(k)) {
        double phiAvg = 0;
        int count = 0;
        for (int idx = 0; idx < TStencil::Q; idx++) {
            if (!m_Geometry.isSolid(mv_Neighbors[k * TStencil::Q+idx])) {
                phiAvg += m_OrderParameter.getParameter(mv_Neighbors[k * TStencil::Q+idx]);
                count++;
            }
        }
        phiAvg /= count;
        m_OrderParameter.getParameter(k) = phiAvg - m_Prefactor * (pow(phiAvg,2) - 1.0);
    }
}


template<typename TLattice, typename TStencil>
inline void CubicWetting<TLattice,TStencil>::setTheta(const double theta){
    m_Prefactor = cos(theta) / sqrt(2.0);
}


template<typename TLattice, typename TStencil>
inline void CubicWetting<TLattice,TStencil>::setThetaDegrees(const double theta){
    setTheta(theta / 180.0 * M_PI);
}


template<typename TLattice, typename TStencil>
inline void CubicWetting<TLattice,TStencil>::communicatePostProcess(){
    m_Data.communicate(m_OrderParameter);
}
