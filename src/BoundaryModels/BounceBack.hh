#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>



template<typename lattice>
class BounceBack : public BoundaryBase {
    public:

        inline void compute(auto& m_Distribution, const int k, const int idx) const;

    private:

};

template<typename lattice>
inline void BounceBack<lattice>::compute(auto& m_Distribution, const int k, const int idx) const {

    m_Distribution.getDistributionPointer(m_Distribution.streamIndex(k,idx))[idx] = m_Distribution.getDistributionPointer(k)[m_Distribution.getOpposite(idx)];

}