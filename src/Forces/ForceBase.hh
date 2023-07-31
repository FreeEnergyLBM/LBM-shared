#pragma once
#include<iostream>
#include "../Forcing.hh"

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

//template<template class objtype, class TMethod,size_t... Is>
//auto constexpr ForceStencilTupleType(std::index_sequence<Is...> is){
//    std::tuple<objtype<TMethod::ForceStencils()[Is]>...> t = std::tie(objtype<TMethod::ForceStencils()[Is]>::getInstance()...);
//    return t;
//}

template<class TMethod = Guo>
class ForceBase{
    
    public:

        using Method = TMethod;

        template<class TTraits>
        inline double computeXYZ(int xyz, int k); //Return force at lattice point k in direction xyz

        template<class TTraits>
        inline double computeQ(int Q, int k); //Return force at lattice point k in direction xyz

        template<class TTraits>
        inline void precompute(int k); //Perform any neccessary computations before force is computed

        template<class TTraits>
        inline double computeDensitySource(int k); //Calculate any possible source/correction term for density

        template<class TTraits>
        inline double computeVelocitySource(int xyz, int k); //Calculate any possible source/correction term for
                                                           //velocity

        template<class TTraits>
        inline void postprocess(int k); //Perform any necessary postprocessing

        template<class TTraits>
        inline void communicatePostProcess(); //Perform any necessary postprocessing

        template<class TTraits>
        inline void communicatePrecompute(); //Perform any necessary postprocessing

    private:

};

template<class TMethod>
template<class TTraits>
inline double ForceBase<TMethod>::computeXYZ(int xyz, int k) {

    return 0;

}

template<class TMethod>
template<class TTraits>
inline double ForceBase<TMethod>::computeQ(int Q, int k) {

    return 0;

}

template<class TMethod>
template<class TTraits>
inline void ForceBase<TMethod>::precompute(int k) { //Not necessary
    
}

template<class TMethod>
template<class TTraits>
inline void ForceBase<TMethod>::postprocess(int k) { //Not necessary
    
}

template<class TMethod>
template<class TTraits>
inline double ForceBase<TMethod>::computeDensitySource(int k) { //Not necessary

    return 0.0;

}

template<class TMethod>
template<class TTraits>
inline double ForceBase<TMethod>::computeVelocitySource(int xyz, int k) { //Need to correct velocity

    return 0;
    
}

template<class TMethod>
template<class TTraits>
inline void ForceBase<TMethod>::communicatePrecompute() { //Not necessary
    
}

template<class TMethod>
template<class TTraits>
inline void ForceBase<TMethod>::communicatePostProcess() { //Not necessary
    
}