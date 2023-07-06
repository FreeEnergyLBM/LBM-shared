#pragma once
#include<iostream>
#include "../Forcing.hh"

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

//template<template class objtype, class T_method,size_t... Is>
//auto constexpr ForceStencilTupleType(std::index_sequence<Is...> is){
//    std::tuple<objtype<T_method::ForceStencils()[Is]>...> t = std::tie(objtype<T_method::ForceStencils()[Is]>::getInstance()...);
//    return t;
//}

template<class T_method = Guo>
class ForceBase{
    
    public:

        using Method = T_method;

        template<class T_traits>
        inline double computeXYZ(int xyz, int k); //Return force at lattice point k in direction xyz

        template<class T_traits>
        inline double computeQ(int Q, int k); //Return force at lattice point k in direction xyz

        template<class T_traits>
        inline void precompute(int k); //Perform any neccessary computations before force is computed

        template<class T_traits>
        inline double computeDensitySource(int k); //Calculate any possible source/correction term for density

        template<class T_traits>
        inline double computeVelocitySource(int xyz, int k); //Calculate any possible source/correction term for
                                                           //velocity

        template<class T_traits>
        inline void postprocess(int k); //Perform any necessary postprocessing

        template<class T_traits>
        inline void communicatePostProcess(); //Perform any necessary postprocessing

        template<class T_traits>
        inline void communicatePrecompute(); //Perform any necessary postprocessing

    private:

};

template<class T_method>
template<class T_traits>
inline double ForceBase<T_method>::computeXYZ(int xyz, int k) {

    return 0;

}

template<class T_method>
template<class T_traits>
inline double ForceBase<T_method>::computeQ(int Q, int k) {

    return 0;

}

template<class T_method>
template<class T_traits>
inline void ForceBase<T_method>::precompute(int k) { //Not necessary
    
}

template<class T_method>
template<class T_traits>
inline void ForceBase<T_method>::postprocess(int k) { //Not necessary
    
}

template<class T_method>
template<class T_traits>
inline double ForceBase<T_method>::computeDensitySource(int k) { //Not necessary

    return 0.0;

}

template<class T_method>
template<class T_traits>
inline double ForceBase<T_method>::computeVelocitySource(int xyz, int k) { //Need to correct velocity

    return 0;
    
}

template<class T_method>
template<class T_traits>
inline void ForceBase<T_method>::communicatePrecompute() { //Not necessary
    
}

template<class T_method>
template<class T_traits>
inline void ForceBase<T_method>::communicatePostProcess() { //Not necessary
    
}