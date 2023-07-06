#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class BounceBack : public BoundaryBase {
    public:

        template<class T_traits, class T_distributiontype>
        inline void compute(T_distributiontype& m_Distribution, int k, int idx);

    private:

};

template<class T_traits, class T_distributiontype>
inline void BounceBack::compute(T_distributiontype& distribution, int k, int idx) {

    distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = distribution.getDistributionPointer(k)[distribution.getOpposite(idx)];

}