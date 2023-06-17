#pragma once

template<class lattice, class stencil>
struct CentralXYZ{

    Data_Base<lattice, stencil> m_Data;

    template<class parameter>
    inline double computeFirstDerivative(const parameter& val, const int direciton, const int k);

    template<class parameter>
    inline double computeLaplacian(const parameter& val, const int k);

    Geometry<lattice> m_Geometry;
    
};

template<class lattice, class stencil>
template<class parameter>
inline double CentralXYZ<lattice, stencil>::computeFirstDerivative(const parameter& val, const int direction, const int k) {

    double gradientsum=0;

    for (int idx = 0; idx <stencil::Q; idx++) {

            gradientsum += 1.0 / stencil::Cs2 * stencil::Weights[idx] * stencil::Ci_xyz(direction)[idx] * 0.5 * (val.getParameter(m_Data.getNeighbors()[k * stencil::Q+idx]));
        
    }

    return gradientsum;

}

template<class lattice, class stencil>
template<class parameter>
inline double CentralXYZ<lattice, stencil>::computeLaplacian(const parameter& val, const int k) {

    double laplaciansum=0;

    for (int idx = 1; idx <stencil::Q; idx++) {

            laplaciansum += 1.0 / stencil::Cs2 * stencil::Weights[idx] * 2 * (val.getParameter(m_Data.getNeighbors()[k * stencil::Q + stencil::Opposites[idx]]) - val.getParameter(k));

    }
    return laplaciansum;
}
