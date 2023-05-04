#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

class BodyForce{
    public:

        template< template< class, class > class data, template< class, int > class parallel, int lx, int ly, int lz = 1 >
        BodyForce( LatticeProperties< data, parallel, lx, ly, lz >& properties, const double magnitude = 0 )
            : m_DT( properties.m_DT ),
              m_Magnitude( magnitude ),
              m_Density( properties ){}

        BodyForce( const BodyForce& other )
            : m_DT( other.m_DT ),
              m_Magnitude( other.m_Magnitude ),
              m_Density( other.m_Density ){}

        double compute( const int xyz, const int k ) const; //Return force at lattice point k in direction xyz

        void precompute( const int k ); //Perform any neccessary computations before force is computed

        double computeDensitySource( const int k ) const; //Calculate any possible source/correction term for density

        double computeVelocitySource( const int xyz, const int k ) const; //Calculate any possible source/correction term for
                                                           //velocity

        void setMagnitude( const double magnitude){
            m_Magnitude=magnitude;
        }

        void postprocess( const int k ); //Perform any necessary postprocessing

    private:

        const double& m_DT;

        double m_Magnitude=0.00000001;//1;

        Density m_Density; //Density

};


double BodyForce::compute( const int xyz, const int k ) const{

    return ( xyz == 0 ) * m_Magnitude * m_Density.getParameter( k ); //Force is just density multiplied by magnitude
                                                                 //in given direction

}

void BodyForce::precompute( const int k ){ //Not necessary
    
}

void BodyForce::postprocess( const int k ){ //Not necessary
    
}

double BodyForce::computeDensitySource( const int k ) const{ //Not necessary

    return 0.0;

}

double BodyForce::computeVelocitySource( const int xyz, const int k) const{ //Need to correct velocity

    return +compute( xyz, k ) * m_DT / ( 2.0 * m_Density.getParameter( k ) );
    
}