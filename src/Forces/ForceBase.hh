#pragma once
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<typename placeholder = void>
class ForceBaseTemplate{
    
    public:

        inline virtual double compute(const int xyz, const int k) const; //Return force at lattice point k in direction xyz

        inline virtual void precompute(const int k); //Perform any neccessary computations before force is computed

        inline virtual double computeDensitySource(const int k) const; //Calculate any possible source/correction term for density

        inline virtual double computeVelocitySource(const int xyz, const int k) const; //Calculate any possible source/correction term for
                                                           //velocity

        inline virtual void postprocess(const int k); //Perform any necessary postprocessing

    private:

};

template<typename placeholder>
inline double ForceBaseTemplate<placeholder>::compute(const int xyz, const int k) const {

    return 0;

}

template<typename placeholder>
inline void ForceBaseTemplate<placeholder>::precompute(const int k) { //Not necessary
    
}

template<typename placeholder>
inline void ForceBaseTemplate<placeholder>::postprocess(const int k) { //Not necessary
    
}

template<typename placeholder>
inline double ForceBaseTemplate<placeholder>::computeDensitySource(const int k) const { //Not necessary

    return 0.0;

}

template<typename placeholder>
inline double ForceBaseTemplate<placeholder>::computeVelocitySource(const int xyz, const int k) const { //Need to correct velocity

    return 0;
    
}

typedef ForceBaseTemplate<> ForceBase;