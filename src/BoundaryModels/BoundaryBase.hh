#pragma once

template< typename placeholder = void >
class BoundaryBaseTemplate{
    public:

        //virtual void compute( auto& m_Distribution, int k, int idx ) const;

        virtual void precompute( const int k );

        virtual double computeDensitySource( const int k ) const;

        virtual double computeVelocitySource( const int xyz, const int k ) const;

        virtual void postprocess( const int k );

    private:


};

template<typename placeholder>
void BoundaryBaseTemplate<placeholder>::precompute( const int k ) {
    
}

template<typename placeholder>
void BoundaryBaseTemplate<placeholder>::postprocess( const int k ) {
    
}

template<typename placeholder>
double BoundaryBaseTemplate<placeholder>::computeDensitySource( const int k ) const {

    return 0;

}

template<typename placeholder>
double BoundaryBaseTemplate<placeholder>::computeVelocitySource( const int xyz,const int k ) const {

    return 0;

}

typedef BoundaryBaseTemplate<> BoundaryBase;