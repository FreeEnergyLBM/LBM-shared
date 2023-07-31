#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class BounceBack : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k, int idx);

    private:

};

template<class TTraits, class TDistributionType>
inline void BounceBack::compute(TDistributionType& distribution, int k, int idx) {

    distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = distribution.getDistributionPointer(k)[distribution.getOpposite(idx)];

}