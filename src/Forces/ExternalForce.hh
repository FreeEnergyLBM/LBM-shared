#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "ForceBase.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<typename placeholder=void>
class BodyForceTemplate : public ForceBase {
    public:

        double compute( const int xyz, const int k ) const override; //Return force at lattice point k in direction xyz

        double computeVelocitySource( const int xyz, const int k ) const override; //Calculate any possible source/correction term for
                                                           //velocity

        void setMagnitude( const double magnitude ){
            m_Magnitude=magnitude;
        }

    private:

        double m_Magnitude=0.00000001;//1;

        Density m_Density; //Density<>

};

template<typename placeholder>
double BodyForceTemplate<placeholder>::compute( const int xyz, const int k ) const{

    return ( xyz == 0 ) * m_Magnitude * m_Density.getParameter( k ); //Force is just density multiplied by magnitude
                                                                 //in given direction

}

template<typename placeholder>
double BodyForceTemplate<placeholder>::computeVelocitySource( const int xyz, const int k) const{ //Need to correct velocity

    return +compute( xyz, k ) * GETPROPERTIES().m_DT / ( 2.0 * m_Density.getParameter( k ) );
    
}

typedef BodyForceTemplate<> BodyForce;