#pragma once

template<typename placeholder=void>
class BoundaryBaseTemplate{
    public:

        //virtual void compute( auto& m_Distribution, int k, int idx ) const;

        virtual void precompute( int k );

        virtual double computeDensitySource( int k ) const;

        virtual double computeVelocitySource( int xyz, int k ) const;

        virtual void postprocess( int k );

    private:


};

template<typename placeholder>
void BoundaryBaseTemplate<placeholder>::precompute( int k ){
    
}

template<typename placeholder>
void BoundaryBaseTemplate<placeholder>::postprocess( int k ){
    
}

template<typename placeholder>
double BoundaryBaseTemplate<placeholder>::computeDensitySource( int k ) const{

    return 0;

}

template<typename placeholder>
double BoundaryBaseTemplate<placeholder>::computeVelocitySource( int xyz,int k ) const{

    return 0;

}

typedef BoundaryBaseTemplate<> BoundaryBase;