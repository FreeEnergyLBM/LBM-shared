#pragma once
#include<iostream>

class AddOnBase{
    
    public:

        template<class T_traits>
        inline double compute(int k) const;

        template<class T_traits>
        inline void communicate();

    private:

};

template<class T_traits>
inline double AddOnBase::compute(int k) const {

    return 0;

}

template<class T_traits>
inline void AddOnBase::communicate() {
    
}