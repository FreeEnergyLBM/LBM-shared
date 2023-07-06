#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "../Forcing.hh"
#include "ForceBase.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<class method=Guo>
class BodyForce : public ForceBase<method> {

    public:

        template<class traits>
        inline double computeXYZ(const int xyz, const int k) const; //Return force at lattice point k in direction xyz

        template<class traits>
        inline double computeQ(const int idx, const int k) const;
        
        template<class traits>
        inline double computeVelocitySource(const int xyz, const int k) const; //Calculate any possible source/correction term for velocity

        inline void setMagnitudeX(const double magnitude);
        inline void setMagnitudeY(const double magnitude);
        inline void setMagnitudeZ(const double magnitude);


    private:

        double m_MagnitudeX = 0;
        double m_MagnitudeY = 0;
        double m_MagnitudeZ = 0;

};

template<class method>
inline void BodyForce<method>::setMagnitudeX(const double magnitude) {

    m_MagnitudeX=magnitude;

}

template<class method>
inline void BodyForce<method>::setMagnitudeY(const double magnitude) {

    m_MagnitudeY=magnitude;

}


template<class method>
inline void BodyForce<method>::setMagnitudeZ(const double magnitude) {

    m_MagnitudeZ=magnitude;

}

template<class method>
template<class traits>
inline double BodyForce<method>::computeXYZ(const int xyz, const int k) const {
    
    if (xyz==0) return m_MagnitudeX * Density<>::get<typename traits::Lattice>(k);
    if constexpr (traits::Lattice::m_NDIM==2){
        return m_MagnitudeY * Density<>::get<typename traits::Lattice>(k);
    }
    else if constexpr (traits::Lattice::m_NDIM==3){
        if(xyz==1) return m_MagnitudeY * Density<>::get<typename traits::Lattice>(k);
        return m_MagnitudeZ * Density<>::get<typename traits::Lattice>(k);
    }
    return 0;
                                                                 //in given direction

}

template<class method>
template<class traits>
inline double BodyForce<method>::computeQ(const int idx, const int k) const {

    if constexpr (traits::Lattice::m_NDIM==2){
        return Density<>::get<typename traits::Lattice>(k) * ( m_MagnitudeX *  traits::Stencil::Ci_xyz(0)[idx]
           + m_MagnitudeY * Density<>::get<typename traits::Lattice>(k) * traits::Stencil::Ci_xyz(1)[idx] );
    }
    else if constexpr (traits::Lattice::m_NDIM==3){
        return Density<>::get<typename traits::Lattice>(k) * ( m_MagnitudeX *  traits::Stencil::Ci_xyz(0)[idx]
           + m_MagnitudeY * Density<>::get<typename traits::Lattice>(k) * traits::Stencil::Ci_xyz(1)[idx]
           + m_MagnitudeZ * Density<>::get<typename traits::Lattice>(k) * traits::Stencil::Ci_xyz(2)[idx] );
    }

}

template<class method>
template<class traits>
inline double BodyForce<method>::computeVelocitySource(const int xyz, const int k) const { //Need to correct velocity

    return +computeXYZ<traits>(xyz, k) * traits::Lattice::m_DT / (2.0);
    
}