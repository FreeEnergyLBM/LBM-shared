#pragma once
#include<iostream>

class AddOnBase{
    
    public:

        template<class TTraits>
        inline void compute(int k);

        template<class TTraits>
        inline void communicate();

    private:

};

template<class TTraits>
inline void AddOnBase::compute(int k){

}

template<class TTraits>
inline void AddOnBase::communicate() {
    
}