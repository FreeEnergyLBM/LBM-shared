#pragma once

template<class lattice, class stencil>
struct CentralXYZNoSolid{

    Data_Base<lattice, stencil, typename std::remove_reference<lattice>::type::template ParallelType<stencil>> m_Data;

    template<class parameter>
    inline double computeFirstDerivative(const parameter& val, const int direciton, const int k);

    template<class parameter>
    inline double computeLaplacian(const parameter& val, const int k);

    Geometry<lattice> m_Geometry;
    
};

template<class lattice, class stencil>
template<class parameter>
inline double CentralXYZNoSolid<lattice, stencil>::computeFirstDerivative(const parameter& val, const int direction, const int k) {

    double gradientsum=0;

    for (int idx = 0; idx <stencil::Q; idx++) {
        
        if ((m_Geometry.isSolid(m_Data.getNeighbors()[k * stencil::Q + idx]))) {

            gradientsum += 1.0 / stencil::Cs2 * stencil::Weights[idx] * stencil::Ci_xyz(direction)[idx] * 0.5 * (val.getParameter(m_Data.getNeighbors()[k * stencil::Q + stencil::Opposites[idx]]));

        }
        else {

            gradientsum += 1.0 / stencil::Cs2 * stencil::Weights[idx] * stencil::Ci_xyz(direction)[idx] * 0.5 * (val.getParameter(m_Data.getNeighbors()[k * stencil::Q+idx]));

        }
        
    }

    return gradientsum;

}

template<class lattice, class stencil>
template<class parameter>
inline double CentralXYZNoSolid<lattice, stencil>::computeLaplacian(const parameter& val, const int k) {

    double laplaciansum=0;

    for (int idx = 1; idx <stencil::Q; idx++) {
    
        if((!m_Geometry.isSolid(m_Data.getNeighbors()[k * stencil::Q + idx]))) {

            laplaciansum += 1.0 / stencil::Cs2 * stencil::Weights[idx] * 2 * (val.getParameter(m_Data.getNeighbors()[k * stencil::Q +idx]) - val.getParameter(k));

        }
        else {

            laplaciansum += 1.0 / stencil::Cs2 * stencil::Weights[idx] * 2 * (val.getParameter(m_Data.getNeighbors()[k * stencil::Q + stencil::Opposites[idx]]) - val.getParameter(k));

        }

    }
    return laplaciansum;
}
