#pragma once
#include <utility>
#include <memory>
#include "../BoundaryModels/BoundaryBase.hh"
#include "../Forces/ForceBase.hh"

template<class lattice=void>
struct DefaultTrait{
    
    using Stencil = std::conditional_t<std::remove_reference<lattice>::type::m_NDIM == 2, D2Q9, D3Q19>; //Here, D refers to the number of cartesian dimensions

    using Boundaries = std::tuple<>;

    using Forces = std::tuple<>;

};

template<class lattice, class traits = DefaultTrait<lattice>>
class ModelBase{ //Inherit from base class to avoid repetition of common
                                                      //calculations
    static_assert(CheckBase<ForceBase, typename traits::Forces>::value, "ERROR: At least one boundary condition chosen is not a boundary class.");
    static_assert(CheckBase<BoundaryBase, typename traits::Boundaries>::value, "ERROR: At least one force chosen is not a force class.");
    public:

        ModelBase()
            : m_Data(),
              m_Distribution(m_Data.getDistributionObject())
        {}

        ModelBase(ModelBase<lattice,traits>& other)
            : m_Data(other.m_Data),
              m_Distribution(other.m_Distribution)
        {}

        inline virtual void precompute(); //Perform any necessary computations before collision

        inline virtual void collide() = 0; //Collision step

        inline virtual void boundaries(); //Boundary calculation

        inline virtual void initialise() = 0; //Initialisation step

        inline virtual void computeMomenta() = 0; //Momenta (density, velocity) calculation

        inline double computeForces(int xyz, int k) const; //Calculate other forces in direction xyz

        inline const std::vector<double>& getDistribution() const; //Return vector of distribution

        template<template<class> class force, int inst = 0>
        inline force<lattice>& getForce() {

            auto forces = get_type<force<lattice>>(mt_Forces);

            return std::get<inst>(forces);

        }

        template<template<class> class boundary, int inst = 0>
        inline boundary<lattice>& getBoundary() {

            auto boundaries = get_type<boundary<lattice>>(mt_Boundaries);

            return std::get<inst>(boundaries);

        }

        typename std::remove_reference<lattice>::type::template DataType<typename traits::Stencil> m_Data; //MOVE THIS TO BASE
        typename std::remove_reference<lattice>::type::template DataType<typename traits::Stencil>::DistributionData& m_Distribution = m_Data.getDistributionObject();
            //Distributions

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions

        typename traits::Forces mt_Forces; //MOVE THIS TO BASE
        typename traits::Boundaries mt_Boundaries; //MOVE THIS TO BASE
        Geometry<lattice> m_Geometry; //MOVE THIS TO BASE
        
        std::vector<double>& distribution = m_Distribution.getDistribution(); //Reference to vector of distributions
        
};

template<class lattice, class traits>
inline const std::vector<double>& ModelBase<lattice,traits>::getDistribution() const {

    return distribution; //Return reference to distribution vector

}

template<class lattice, class traits>
inline void ModelBase<lattice,traits>::precompute() {

    #pragma omp for schedule(guided)
    for (int k = lattice::m_HaloSize; k <lattice::m_N - lattice::m_HaloSize; k++) { //loop over k

        if constexpr(std::tuple_size<typename traits::Forces>::value != 0){ //Check if there is at least one element
                                                                          //in F

            std::apply([k](auto&... forces){//See Algorithm.hh for explanation of std::apply

                (forces.precompute(k),...);

            }, mt_Forces);

        }
        
    }

    #pragma omp master
    {
    m_Distribution.getDistribution().swap(m_Distribution.getDistributionOld()); //swap old and new distributions
                                                                                //before collision
    }
    
}

template<class lattice, class traits>
inline double ModelBase<lattice,traits>::computeForces(int xyz, int k) const {

    if constexpr(std::tuple_size<typename traits::Forces>::value != 0){

        return std::apply([xyz, k](auto&... forces){
                return (forces.compute(xyz, k) + ...);
            }, mt_Forces);

    }
    else return 0;

}

template<class lattice, class traits>
inline void ModelBase<lattice,traits>::boundaries() {

    #pragma omp for schedule(guided)
    for (int k = 0; k <lattice::m_N; k++) { //loop over k

        if constexpr(std::tuple_size<typename traits::Boundaries>::value != 0) { //Check if there are any boundary
                                                                              //models

            std::apply([this, k](auto&... boundaries) {
                // Make this a sub function for readability
                for (int idx = 0; idx <traits::Stencil::Q; idx++) {

                    if(m_Geometry.isSolid(k) && !m_Geometry.isSolid(m_Distribution.streamIndex(k, idx))) {

                        (boundaries.compute(this -> m_Distribution, k, idx) , ...);

                    }

                }

            }, mt_Boundaries);
            
        }

    }
    
}
