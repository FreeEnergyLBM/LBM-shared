#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template< class gradientstencil >
class OrderParameterGradients{
    public:
        template< template< class, class > class data, template< class, int > class parallel, int lx, int ly, int lz = 1 >
        OrderParameterGradients( LatticeProperties< data, parallel, lx, ly, lz >& properties ) 
         : m_GradientStencil( properties ),
           m_GradOrderParameter( properties ),
           m_LaplacianOrderParameter( properties ),
           m_OrderParameter( properties ),
           m_NDIM( properties.m_NDIM )
        {}

        OrderParameterGradients( const OrderParameterGradients<gradientstencil>& other ): m_GradientStencil( other.m_GradientStencil ),
                                                                                        m_GradOrderParameter( other.m_GradOrderParameter ),
                                                                                        m_LaplacianOrderParameter( other.m_LaplacianOrderParameter ),
                                                                                        m_OrderParameter( other.m_OrderParameter ),
                                                                                        NDIM( other.NDIM ){}

        double compute( int xyz, int k ) const; //Return force at lattice point k in direction xyz

        void precompute( int k ); //Perform any neccessary computations before force is computed

        double computeDensitySource( int k ) const; //Calculate any possible source/correction term for density

        double computeVelocitySource( int xyz,int k ) const; //Calculate any possible source/correction term for
                                                           //velocity

        void postprocess( int k ); //Perform any necessary postprocessing

        gradientstencil m_GradientStencil;

        GradientOrderParameter m_GradOrderParameter;

        LaplacianOrderParameter m_LaplacianOrderParameter;

        OrderParameter m_OrderParameter;

        const int& m_NDIM;

};

template< class gradientstencil >
double OrderParameterGradients< gradientstencil >::compute( int xyz,int k ) const{

    return 0;

}

template< class gradientstencil >
void OrderParameterGradients< gradientstencil >::precompute( int k ){ //Not necessary
    
    for( int xyz = 0; xyz < m_NDIM; xyz++ ) m_GradOrderParameter.getParameterPointer( k )[ xyz ] = m_GradientStencil.computeFirstDerivative( m_OrderParameter, xyz, k );
    m_LaplacianOrderParameter.getParameter( k ) = m_GradientStencil.computeLaplacian( m_OrderParameter, k );

}

template< class gradientstencil >
void OrderParameterGradients< gradientstencil >::postprocess( int k ){ //Not necessary
    
}

template< class gradientstencil >
double OrderParameterGradients< gradientstencil >::computeDensitySource( int k ) const{ //Not necessary

    return 0.0;

}

template< class gradientstencil >
double OrderParameterGradients< gradientstencil >::computeVelocitySource( int xyz,int k ) const{ //Need to correct velocity

    return 0;
    
}
