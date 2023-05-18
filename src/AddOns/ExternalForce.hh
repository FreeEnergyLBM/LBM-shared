#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "AddOnBase.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<class lattice>
class BodyForce : public AddOnBase {

    public:

        inline double compute(const int xyz, const int k) const override; //Return force at lattice point k in direction xyz

        inline double computeVelocitySource(const int xyz, const int k) const override; //Calculate any possible source/correction term for velocity

        inline void setMagnitudeX(const double magnitude);
        inline void setMagnitudeY(const double magnitude);
        inline void setMagnitudeZ(const double magnitude);


    private:

        double m_MagnitudeX = 0;
        double m_MagnitudeY = 0;
        double m_MagnitudeZ = 0;

        Density<lattice> m_Density; //Density<>

};

template<typename lattice>
inline void BodyForce<lattice>::setMagnitudeX(const double magnitude) {

    m_MagnitudeX=magnitude;

}

template<typename lattice>
inline void BodyForce<lattice>::setMagnitudeY(const double magnitude) {

    m_MagnitudeY=magnitude;

}


template<typename lattice>
inline void BodyForce<lattice>::setMagnitudeZ(const double magnitude) {

    m_MagnitudeZ=magnitude;

}

template<typename lattice>
inline double BodyForce<lattice>::compute(const int xyz, const int k) const {

    if (xyz==0) return m_MagnitudeX * m_Density.getParameter(k);
    if constexpr (lattice::m_NDIM==2){
        return m_MagnitudeY * m_Density.getParameter(k);
    }
    else if constexpr (lattice::m_NDIM==3){
        if(xyz==1) return m_MagnitudeY * m_Density.getParameter(k);
        return m_MagnitudeZ * m_Density.getParameter(k);
    }
    return 0;
    
                                                                 //in given direction

}

template<class lattice>
inline double BodyForce<lattice>::computeVelocitySource(const int xyz, const int k) const { //Need to correct velocity

    return +compute(xyz, k) * lattice::m_DT / (2.0 * m_Density.getParameter(k));
    
}