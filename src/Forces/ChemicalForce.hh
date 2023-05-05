#pragma once
#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "ForceBase.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template< typename placeholder = void >
class ChemicalForceTemplate : public ForceBase {
    
    public:

        double compute( const int xyz, const int k ) const override; //Return force at lattice point k in direction xyz

        void precompute( const int k ) override; //Perform any neccessary computations before force is computed

        double computeVelocitySource( const int xyz,const int k ) const override; //Calculate any possible source/correction term for
                                                           //velocity

    private:

        double m_A=0.00015;

        double m_Kappa=0.0003;

        ChemicalPotential m_ChemicalPotential;

        GradientOrderParameter m_GradOrderParameter;

        LaplacianOrderParameter m_LaplacianOrderParameter;

        OrderParameter m_OrderParameter;

        Density m_Density;

};

template< typename placeholder >
double ChemicalForceTemplate< placeholder >::compute( const int xyz, const int k ) const {

    return m_ChemicalPotential.getParameter( k ) * m_GradOrderParameter.getParameterPointer( k )[ xyz ];

}

template< typename placeholder >
void ChemicalForceTemplate< placeholder >::precompute( const int k ){ //Not necessary

    double orderparamcubed=m_OrderParameter.getParameter( k ) * m_OrderParameter.getParameter( k ) * m_OrderParameter.getParameter( k );
    m_ChemicalPotential.getParameter( k ) = -m_A * m_OrderParameter.getParameter( k ) + m_A * orderparamcubed - m_Kappa * m_LaplacianOrderParameter.getParameter( k );

}

template< typename placeholder >
double ChemicalForceTemplate< placeholder >::computeVelocitySource( const int xyz, const int k ) const{ //Need to correct velocity

    return +compute(xyz,k) * GETPROPERTIES().m_DT / ( 2.0 * m_Density.getParameter( k ) );
    
}

typedef ChemicalForceTemplate<> ChemicalForce;