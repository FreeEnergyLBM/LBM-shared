#pragma once
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

class AddOnBase{
    
    public:

        template<class traits>
        inline double compute(const int k) const;

        template<class traits>
        inline void communicate();

    private:

};

template<class traits>
inline double AddOnBase::compute(const int k) const {

    return 0;

}

template<class traits>
inline void AddOnBase::communicate() { //Not necessary
    
}