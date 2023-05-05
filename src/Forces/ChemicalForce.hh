#pragma once
#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<typename placeholder=void>
class ChemicalForceTemplate{
    public:

        double compute( int xyz, int k ) const; //Return force at lattice point k in direction xyz

        void precompute( int k ); //Perform any neccessary computations before force is computed

        double computeDensitySource( int k ) const; //Calculate any possible source/correction term for density

        double computeVelocitySource( int xyz,int k ) const; //Calculate any possible source/correction term for
                                                           //velocity

        void postprocess( int k ); //Perform any necessary postprocessing

    private:

        double m_A=0.00015;

        double m_Kappa=0.0003;

        ChemicalPotential m_ChemicalPotential;

        GradientOrderParameter m_GradOrderParameter;

        LaplacianOrderParameter m_LaplacianOrderParameter;

        OrderParameter m_OrderParameter;

        Density m_Density; //Density

};

template<typename placeholder>
double ChemicalForceTemplate<placeholder>::compute( int xyz, int k ) const{

    return m_ChemicalPotential.getParameter( k ) * m_GradOrderParameter.getParameterPointer( k )[ xyz ];

}

template<typename placeholder>
void ChemicalForceTemplate<placeholder>::precompute( int k ){ //Not necessary

    m_ChemicalPotential.getParameter( k ) = -m_A * m_OrderParameter.getParameter( k ) + m_A * m_OrderParameter.getParameter( k ) * m_OrderParameter.getParameter( k ) * m_OrderParameter.getParameter( k ) - m_Kappa * m_LaplacianOrderParameter.getParameter( k );

}

template<typename placeholder>
void ChemicalForceTemplate<placeholder>::postprocess( int k ){ //Not necessary
    
}

template<typename placeholder>
double ChemicalForceTemplate<placeholder>::computeDensitySource( int k ) const{ //Not necessary

    return 0.0;

}

template<typename placeholder>
double ChemicalForceTemplate<placeholder>::computeVelocitySource( int xyz, int k ) const{ //Need to correct velocity

    return +compute(xyz,k) * GETPROPERTIES().m_DT / ( 2.0 * m_Density.getParameter( k ) );
    
}

typedef ChemicalForceTemplate<> ChemicalForce;