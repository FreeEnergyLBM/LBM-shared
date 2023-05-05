#pragma once
#include <utility>

template< class properties >
struct DefaultTrait{
    
    using Stencil = std::conditional_t< std::remove_reference< properties >::type::m_NDIM == 2, D2Q9, D3Q19 >; //Here, D refers to the number of cartesian dimensions

    using Boundaries = std::tuple<>;

    using Forces = std::tuple<>;

    using Properties = properties;

};

template<class traits = DefaultTrait< decltype( GETPROPERTIES() ) > >
class ModelBase{ //Inherit from base class to avoid repetition of common
                                                      //calculations
                                                      
    public:

        ModelBase()
            : m_Data(),
              m_Distribution( m_Data.getDistributionObject() )
        {}

        ModelBase( ModelBase< traits >& other )
            : m_Data( other.m_Data ),
              m_Distribution( other.m_Distribution )
        {}

        virtual void precompute(); //Perform any necessary computations before collision

        virtual void collide() = 0; //Collision step

        virtual void boundaries(); //Boundary calculation

        virtual void initialise() = 0; //Initialisation step

        virtual void computeMomenta() = 0; //Momenta (density, velocity) calculation

        double computeForces( int xyz, int k ) const; //Calculate other forces in direction xyz

        const std::vector< double >& getDistribution() const; //Return vector of distribution

        template< class force, int inst = 0 >
        force& getForce() {

            auto forces = get_type< force >( mt_Forces );

            return std::get< inst >( forces );

        }

        template< class boundary, int inst = 0 >
        boundary& getBoundary() {

            auto boundaries = get_type< boundary >( mt_Boundaries );

            return std::get< inst >( boundaries );

        }

        typename std::remove_reference< typename traits::Properties >::type::template DataType< typename traits::Stencil > m_Data; //MOVE THIS TO BASE
        typename std::remove_reference< typename traits::Properties >::type::template DataType< typename traits::Stencil >::DistributionData& m_Distribution = m_Data.getDistributionObject();
            //Distributions

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions

        typename traits::Forces mt_Forces; //MOVE THIS TO BASE
        typename traits::Boundaries mt_Boundaries; //MOVE THIS TO BASE
        Geometry m_Geometry; //MOVE THIS TO BASE

        vector< double >& distribution = m_Distribution.getDistribution(); //Reference to vector of distributions
        
};

template< class traits >
const std::vector< double >& ModelBase< traits >::getDistribution() const {

    return distribution; //Return reference to distribution vector

}

template< class traits >
void ModelBase< traits >::precompute() {

    #pragma omp parallel for schedule( guided )
    for ( int k = GETPROPERTIES().m_HaloSize; k < GETPROPERTIES().m_N - GETPROPERTIES().m_HaloSize; k++ ) { //loop over k

        if constexpr( std::tuple_size< typename traits::Forces >::value != 0 ){ //Check if there is at least one element
                                                                          //in F

            std::apply([ k ]( auto&... forces ){//See Algorithm.hh for explanation of std::apply

                (forces.precompute( k ),...);

            }, mt_Forces );

        }
        
    }

    m_Distribution.getDistribution().swap( m_Distribution.getDistributionOld() ); //swap old and new distributions
                                                                                //before collision

}

template< class traits >
double ModelBase< traits >::computeForces( int xyz, int k ) const {

    if constexpr( std::tuple_size< typename traits::Forces >::value != 0 ){

        return std::apply([ xyz, k ]( auto&... forces ){
                return ( forces.compute( xyz, k ) + ... );
            }, mt_Forces );

    }
    else return 0;

}

template< class traits >
void ModelBase< traits >::boundaries() {

    #pragma omp parallel for schedule( guided )
    for ( int k = 0; k < GETPROPERTIES().m_N; k++ ) { //loop over k

        if constexpr( std::tuple_size< typename traits::Boundaries >::value != 0 ) { //Check if there are any boundary
                                                                              //models

            std::apply( [ this, k ]( auto&... boundaries ) {
                // Make this a sub function for readability
                for ( int idx = 0; idx < traits::Stencil::Q; idx++ ) {

                    if( m_Geometry.isSolid( k ) && !m_Geometry.isSolid( m_Distribution.streamIndex( k, idx ) ) ) {

                        ( boundaries.compute( this -> m_Distribution, k, idx ) , ... );

                    }

                }

            }, mt_Boundaries );
            
        }

    }
    
}