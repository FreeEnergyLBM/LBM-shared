#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "AddOnBase.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<class lattice, class gradientstencil>
class OrderParameterGradients : public AddOnBase {

    public:

        inline void precompute(const int k) override; //Perform any neccessary computations before force is computed

    private:

        gradientstencil m_GradientStencil;

        GradientOrderParameter<lattice> m_GradOrderParameter;

        LaplacianOrderParameter<lattice> m_LaplacianOrderParameter;

        OrderParameter<lattice> m_OrderParameter;

};


template<class lattice, class gradientstencil>
inline void OrderParameterGradients<lattice, gradientstencil>::precompute(const int k) { //Not necessary
    
    for(int xyz = 0; xyz <lattice::m_NDIM; xyz++) m_GradOrderParameter.getParameterPointer(k)[xyz] = m_GradientStencil.computeFirstDerivative(m_OrderParameter, xyz, k);

    m_LaplacianOrderParameter.getParameter(k) = m_GradientStencil.computeLaplacian(m_OrderParameter, k);

}