#pragma once
#include <iostream>

#include "../Forcing.hh"

class Model;

template <class TMethod = Guo>
class ForceBase {
   public:
    using Method = TMethod;

    template <class TTraits>
    inline double computeXYZ(int xyz, int k);  // Return force at lattice point k in direction xyz

    template <class TTraits>
    inline double computeQ(int Q, int k);  // Return force at lattice point k in direction q

    template <class TTraits>
    inline double compute(int k);

    template <class TTraits>
    inline void runProcessor(int k);  // Perform any neccessary computations before force is computed

    template <class TTraits>
    inline double computeDensitySource(int k);  // Calculate any possible source/correction term for density

    template <class TTraits>
    inline double computeDensitySourceMultiplicative(
        int k);  // Calculate any possible source/correction term for density

    template <class TTraits>
    inline double computeVelocitySource(int xyz, int k);  // Calculate any possible source/correction term for
                                                          // velocity

    template <class TTraits>
    inline void communicate();  // Perform any necessary postprocessing

    inline void initialise(Model* model) { mModel = model; };

   protected:
    Model* mModel;
};

template <class TMethod>
template <class TTraits>
inline double ForceBase<TMethod>::computeXYZ(int xyz, int k) {
    return 0;
}

template <class TMethod>
template <class TTraits>
inline double ForceBase<TMethod>::computeQ(int Q, int k) {
    return 0;
}

template <class TMethod>
template <class TTraits>
inline double ForceBase<TMethod>::compute(int k) {
    return 0;
}

template <class TMethod>
template <class TTraits>
inline void ForceBase<TMethod>::runProcessor(int k) {  // Not necessary
}

template <class TMethod>
template <class TTraits>
inline double ForceBase<TMethod>::computeDensitySource(int k) {  // Not necessary

    return 0.0;
}

template <class TMethod>
template <class TTraits>
inline double ForceBase<TMethod>::computeDensitySourceMultiplicative(int k) {  // Not necessary

    return 1.0;
}

template <class TMethod>
template <class TTraits>
inline double ForceBase<TMethod>::computeVelocitySource(int xyz, int k) {  // Need to correct velocity

    return 0;
}

template <class TMethod>
template <class TTraits>
inline void ForceBase<TMethod>::communicate() {  // Not necessary
}
