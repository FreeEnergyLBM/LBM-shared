#pragma once
#include "../Parameters.hh"
#include "AddOnBase.hh"
#include<iostream>
#include<math.h>



template<typename lattice, typename stencil>
class LinearWetting : public AddOnBase {
    public:

        inline void postprocess(const int k);

        inline void setTheta(const double theta);

        inline void setThetaDegrees(const double theta);

        inline void setPrefactor(const double prefactor);

        inline void setPrefactor(const double A, const double kappa);

        inline void communicatePostProcess();

    private:

        inline double calcOmega(const double theta);

        OrderParameter<lattice> m_OrderParameter;
        Geometry<lattice> m_Geometry;

        Data_Base<lattice, stencil, typename lattice::template ParallelType<stencil>> m_Data;

        double m_Theta = M_PI / 2.0;
        double m_Omega;
        double m_Prefactor = 0.0;
        const std::vector<int> mv_Neighbors = m_Data.getNeighbors();

};

template<typename lattice, typename stencil>
inline void LinearWetting<lattice,stencil>::postprocess(const int k) {

    if (m_Geometry.isSolid(k)) {
        
        double wettingsum = 0;
        int count = 0;

        for (int idx = 0; idx < stencil::Q; idx++) {
                
            if (!m_Geometry.isSolid(mv_Neighbors[k * stencil::Q+idx])) {

                wettingsum += m_Prefactor * m_Omega + m_OrderParameter.getParameter(mv_Neighbors[k * stencil::Q+idx]);
                count++;

            }

        }

        m_OrderParameter.getParameter(k) = wettingsum/((double)count);

    }
    
}

template<typename lattice, typename stencil>
inline void LinearWetting<lattice,stencil>::setTheta(const double theta){

    m_Theta = theta;
    m_Omega = calcOmega(m_Theta);

}

template<typename lattice, typename stencil>
inline void LinearWetting<lattice,stencil>::setThetaDegrees(const double theta){

    m_Theta = M_PI * theta / 180.0;
    m_Omega = calcOmega(m_Theta);

}

template<typename lattice, typename stencil>
inline double LinearWetting<lattice,stencil>::calcOmega(const double theta){

    double alpha = acos(sin(theta) * sin(theta));
    return 2 * (((M_PI / 2.0 - theta) >= 0) - ((M_PI / 2.0 - theta) < 0)) * sqrt(cos(alpha / 3.0) * (1.0 - cos(alpha / 3.0)));

}

template<typename lattice, typename stencil>
inline void LinearWetting<lattice,stencil>::setPrefactor(const double prefactor){

    m_Prefactor=prefactor;

}

template<typename lattice, typename stencil>
inline void LinearWetting<lattice,stencil>::setPrefactor(const double A, const double kappa){

    m_Prefactor=sqrt(A/(2.0*kappa));

}

template<typename lattice, typename stencil>
inline void LinearWetting<lattice,stencil>::communicatePostProcess(){

    #ifdef MPIPARALLEL
    #pragma omp master
    {
    m_Data.communicate(m_OrderParameter);
    }
    #endif
    
}