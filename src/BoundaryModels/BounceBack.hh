#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>



template<typename lattice>
class BounceBack : public BoundaryBase {
    public:

        template<class distribution>
        inline void compute(distribution& m_Distribution, const int k, const int idx) const;

    private:

};

template<typename lattice>
template<class distribution>
inline void BounceBack<lattice>::compute(distribution& m_Distribution, const int k, const int idx) const {

    m_Distribution.getDistributionPointer(m_Distribution.streamIndex(k,idx))[idx] = m_Distribution.getDistributionPointer(k)[m_Distribution.getOpposite(idx)];

}