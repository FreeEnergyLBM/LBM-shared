#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "../Forcing.hh"
#include "ForceBase.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<class TMethod = Guo>
class BodyForce : public ForceBase<TMethod> {

    public:

        template<class TTraits>
        inline double computeXYZ(int xyz, int k); //Return force at lattice point k in direction xyz

        template<class TTraits>
        inline double computeQ(int idx, int k);
        
        template<class TTraits>
        inline double computeVelocitySource(int xyz, int k); //Calculate any possible source/correction term for velocity

        inline void setMagnitudeX(double magnitude);
        inline void setMagnitudeY(double magnitude);
        inline void setMagnitudeZ(double magnitude);

        // These functions activate gravity in a desired (y or z) direction.
        // the Gravity term is only applied to the c_1: c_1 = 0.5 * rho * (1 + phi)
        inline void activateGravityY() { gravityY = true; }
        inline void activateGravityZ() { gravityZ = true; }

    private:

        double mMagnitudeX = 0;
        double mMagnitudeY = 0;
        double mMagnitudeZ = 0;

        bool gravityY = false;
        bool gravityZ = false;
};

template<class TMethod>
inline void BodyForce<TMethod>::setMagnitudeX(double magnitude) {

    mMagnitudeX = magnitude;

}

template<class TMethod>
inline void BodyForce<TMethod>::setMagnitudeY(double magnitude) {

    mMagnitudeY = magnitude;

}


template<class TMethod>
inline void BodyForce<TMethod>::setMagnitudeZ(double magnitude) {

    mMagnitudeZ = magnitude;

}

template<class TMethod>
template<class TTraits>
inline double BodyForce<TMethod>::computeXYZ(int xyz, int k) {

    using Lattice = typename TTraits::Lattice;

    // Why is Density<>::get<typename TTraits::Lattice>(k) not used and instead a constant density is used?
    double density = 1.;    
    
    // 2D cases
    if constexpr (Lattice::NDIM == 2){
        if (xyz == 0) return mMagnitudeX * density;
        
        if (gravityY == true)
            return mMagnitudeY * Density<>::get<typename TTraits::Lattice>(k) * (1.0 + OrderParameter<>::get<typename TTraits::Lattice>(k)) / 2.0;
        else
        return mMagnitudeY * density;
    }

    // 3D cases
    else if constexpr (Lattice::NDIM == 3){
        if (xyz == 0) return mMagnitudeX * density;

        if (xyz == 1) {
        if (gravityY == true)
            return mMagnitudeY * Density<>::get<typename TTraits::Lattice>(k) * (1.0 + OrderParameter<>::get<typename TTraits::Lattice>(k)) / 2.0;
        else return mMagnitudeY * density;}
        
        if (gravityZ == true)
            return mMagnitudeZ * Density<>::get<typename TTraits::Lattice>(k) * (1.0 + OrderParameter<>::get<typename TTraits::Lattice>(k)) / 2.0;
        else return mMagnitudeZ * density;
    }

    return mMagnitudeX * density;

}

template<class TMethod>
template<class TTraits>
inline double BodyForce<TMethod>::computeQ(int idx, int k) {

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    //double& density = Density<>::get<Lattice>(k);
    double density = 1.;

    if constexpr (Lattice::NDIM == 2){
        return density * ( mMagnitudeX * Stencil::Ci_xyz(0)[idx]
           + mMagnitudeY *  Stencil::Ci_xyz(1)[idx] );
    }

    else if constexpr (Lattice::NDIM == 3){
        return density * ( mMagnitudeX * Stencil::Ci_xyz(0)[idx]
           + mMagnitudeY * Stencil::Ci_xyz(1)[idx]
           + mMagnitudeZ * Stencil::Ci_xyz(2)[idx] );
    }

    return density * ( mMagnitudeX * Stencil::Ci_xyz(0)[idx] );

}

template<class TMethod>
template<class TTraits>
inline double BodyForce<TMethod>::computeVelocitySource(int xyz, int k) {

    return +computeXYZ<TTraits>(xyz, k) * TTraits::Lattice::DT / (2.0);
    
}