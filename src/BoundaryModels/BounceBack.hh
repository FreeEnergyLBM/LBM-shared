#pragma once
#include "../Parameters.hh"
#include<iostream>

template<typename placeholder=void>
class BounceBackTemplate{
    public:

        void compute( auto& m_Distribution, int k, int idx ) const;

        void precompute( int k );

        double computeDensitySource( int k ) const;

        double computeVelocitySource( int xyz, int k ) const;

        void postprocess( int k );

    private:


};

template<typename placeholder>
void BounceBackTemplate<placeholder>::compute( auto& m_Distribution, int k, int idx) const{
    m_Distribution.getDistributionPointer( m_Distribution.streamIndex(k,idx) )[ idx ] = m_Distribution.getDistributionPointer( k )[ m_Distribution.getOpposite( idx ) ];
}

template<typename placeholder>
void BounceBackTemplate<placeholder>::precompute( int k ){
    
}

template<typename placeholder>
void BounceBackTemplate<placeholder>::postprocess( int k ){
    
}

template<typename placeholder>
double BounceBackTemplate<placeholder>::computeDensitySource( int k ) const{

    return 0;

}

template<typename placeholder>
double BounceBackTemplate<placeholder>::computeVelocitySource( int xyz,int k ) const{

    return 0;

}

typedef BounceBackTemplate<> BounceBack;