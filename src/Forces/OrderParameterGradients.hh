#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "ForceBase.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template< class gradientstencil >
class OrderParameterGradients : public ForceBase {

    public:

        void precompute( const int k ) override; //Perform any neccessary computations before force is computed

    private:

        gradientstencil m_GradientStencil;

        GradientOrderParameter m_GradOrderParameter;

        LaplacianOrderParameter m_LaplacianOrderParameter;

        OrderParameter m_OrderParameter;

};


template< class gradientstencil >
void OrderParameterGradients< gradientstencil >::precompute( const int k ) { //Not necessary
    
    for( int xyz = 0; xyz < GETPROPERTIES().m_NDIM; xyz++ ) m_GradOrderParameter.getParameterPointer( k )[ xyz ] = m_GradientStencil.computeFirstDerivative( m_OrderParameter, xyz, k );

    m_LaplacianOrderParameter.getParameter( k ) = m_GradientStencil.computeLaplacian( m_OrderParameter, k );

}