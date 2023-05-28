#pragma once
#include<iostream>
#include "../Forcing.hh"

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

class AddOnBase{
    
    public:

        using Method = ForcingNone;

        inline virtual double computeXYZ(const int xyz, const int k) const; //Return force at lattice point k in direction xyz

        inline virtual double computeQ(const int Q, const int k) const; //Return force at lattice point k in direction xyz

        inline virtual void precompute(const int k); //Perform any neccessary computations before force is computed

        inline virtual double computeDensitySource(const int k) const; //Calculate any possible source/correction term for density

        inline virtual double computeVelocitySource(const int xyz, const int k) const; //Calculate any possible source/correction term for
                                                           //velocity

        inline virtual void postprocess(const int k); //Perform any necessary postprocessing

        inline virtual void communicatePostProcess(); //Perform any necessary postprocessing

        inline virtual void communicatePrecompute(); //Perform any necessary postprocessing

    private:

};

inline double AddOnBase::computeXYZ(const int xyz, const int k) const {

    return 0;

}

inline double AddOnBase::computeQ(const int Q, const int k) const {

    return 0;

}

inline void AddOnBase::precompute(const int k) { //Not necessary
    
}

inline void AddOnBase::postprocess(const int k) { //Not necessary
    
}

inline double AddOnBase::computeDensitySource(const int k) const { //Not necessary

    return 0.0;

}

inline double AddOnBase::computeVelocitySource(const int xyz, const int k) const { //Need to correct velocity

    return 0;
    
}

inline void AddOnBase::communicatePrecompute() { //Not necessary
    
}

inline void AddOnBase::communicatePostProcess() { //Not necessary
    
}