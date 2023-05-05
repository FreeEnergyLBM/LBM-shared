#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<typename placeholder=void>
class BodyForceTemplate{
    public:

        double compute( const int xyz, const int k ) const; //Return force at lattice point k in direction xyz

        void precompute( const int k ); //Perform any neccessary computations before force is computed

        double computeDensitySource( const int k ) const; //Calculate any possible source/correction term for density

        double computeVelocitySource( const int xyz, const int k ) const; //Calculate any possible source/correction term for
                                                           //velocity

        void setMagnitude( const double magnitude ){
            m_Magnitude=magnitude;
        }

        void postprocess( const int k ); //Perform any necessary postprocessing

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
void BodyForceTemplate<placeholder>::precompute( const int k ){ //Not necessary
    
}

template<typename placeholder>
void BodyForceTemplate<placeholder>::postprocess( const int k ){ //Not necessary
    
}

template<typename placeholder>
double BodyForceTemplate<placeholder>::computeDensitySource( const int k ) const{ //Not necessary

    return 0.0;

}

template<typename placeholder>
double BodyForceTemplate<placeholder>::computeVelocitySource( const int xyz, const int k) const{ //Need to correct velocity

    return +compute( xyz, k ) * GETPROPERTIES().m_DT / ( 2.0 * m_Density.getParameter( k ) );
    
}

typedef BodyForceTemplate<> BodyForce;