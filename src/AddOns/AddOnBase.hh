#pragma once
#include<iostream>

class AddOnBase{
    
    public:

        template<class T_traits>
        inline void compute(int k);

        template<class T_traits>
        inline void communicate();

    private:

};

template<class T_traits>
inline void AddOnBase::compute(int k){

}

template<class T_traits>
inline void AddOnBase::communicate() {
    
}