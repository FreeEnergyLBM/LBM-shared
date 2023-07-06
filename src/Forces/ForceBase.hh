#pragma once
#include<iostream>
#include "../Forcing.hh"

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

//template<template class objtype, class method,size_t... Is>
//auto constexpr ForceStencilTupleType(std::index_sequence<Is...> is){
//    std::tuple<objtype<method::ForceStencils()[Is]>...> t = std::tie(objtype<method::ForceStencils()[Is]>::getInstance()...);
//    return t;
//}

template<class method=Guo>
class ForceBase{
    
    public:

        using Method = method;

        template<class traits>
        inline double computeXYZ(const int xyz, const int k) const; //Return force at lattice point k in direction xyz

        template<class traits>
        inline double computeQ(const int Q, const int k) const; //Return force at lattice point k in direction xyz

        template<class traits>
        inline void precompute(const int k); //Perform any neccessary computations before force is computed

        template<class traits>
        inline double computeDensitySource(const int k) const; //Calculate any possible source/correction term for density

        template<class traits>
        inline double computeVelocitySource(const int xyz, const int k) const; //Calculate any possible source/correction term for
                                                           //velocity

        template<class traits>
        inline void postprocess(const int k); //Perform any necessary postprocessing

        template<class traits>
        inline void communicatePostProcess(); //Perform any necessary postprocessing

        template<class traits>
        inline void communicatePrecompute(); //Perform any necessary postprocessing

    private:

};

template<class method>
template<class traits>
inline double ForceBase<method>::computeXYZ(const int xyz, const int k) const {

    return 0;

}

template<class method>
template<class traits>
inline double ForceBase<method>::computeQ(const int Q, const int k) const {

    return 0;

}

template<class method>
template<class traits>
inline void ForceBase<method>::precompute(const int k) { //Not necessary
    
}

template<class method>
template<class traits>
inline void ForceBase<method>::postprocess(const int k) { //Not necessary
    
}

template<class method>
template<class traits>
inline double ForceBase<method>::computeDensitySource(const int k) const { //Not necessary

    return 0.0;

}

template<class method>
template<class traits>
inline double ForceBase<method>::computeVelocitySource(const int xyz, const int k) const { //Need to correct velocity

    return 0;
    
}

template<class method>
template<class traits>
inline void ForceBase<method>::communicatePrecompute() { //Not necessary
    
}

template<class method>
template<class traits>
inline void ForceBase<method>::communicatePostProcess() { //Not necessary
    
}