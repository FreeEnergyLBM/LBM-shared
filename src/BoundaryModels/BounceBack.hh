#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>

template<typename placeholder = void>
class BounceBackTemplate : public BoundaryBase {
    public:

        void compute(auto& m_Distribution, const int k, const int idx) const;

    private:


};

template<typename placeholder>
void BounceBackTemplate<placeholder>::compute(auto& m_Distribution, const int k, const int idx) const {

    m_Distribution.getDistributionPointer(m_Distribution.streamIndex(k,idx))[idx] = m_Distribution.getDistributionPointer(k)[m_Distribution.getOpposite(idx)];

}

typedef BounceBackTemplate<> BounceBack;