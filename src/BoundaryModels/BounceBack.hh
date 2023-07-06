#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class BounceBack : public BoundaryBase {
    public:

        template<class traits, class distributiontype>
        inline void compute(distributiontype& m_Distribution, const int k, const int idx) const;

    private:

};

template<class traits, class distributiontype>
inline void BounceBack::compute(distributiontype& distribution, const int k, const int idx) const {

    distribution.getDistributionPointer(distribution.streamIndex(k,idx))[idx] = distribution.getDistributionPointer(k)[distribution.getOpposite(idx)];

}