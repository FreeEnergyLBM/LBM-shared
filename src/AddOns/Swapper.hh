#pragma once
#include<iostream>

template<class TParam1, class TParam2, int numdir = 1>
class Swapper : public AddOnBase {
    
    public:

        template<class TTraits>
        inline void compute(int k);

    private:

};

template<class TParam1, class TParam2, int numdir>
template<class TTraits>
inline void Swapper<TParam1, TParam2, numdir>::compute(int k){

    TParam1::template get<typename TTraits::Lattice, numdir>().swap(TParam2::template get<typename TTraits::Lattice, numdir>());

}