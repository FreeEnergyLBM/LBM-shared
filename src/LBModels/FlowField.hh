#pragma once
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include "../BoundaryModels/Boundaries.hh"
#include "../AddOns/AddOns.hh"
#include "../Forces/Forces.hh"
#include "../GradientStencils/GradientStencils.hh"
#include "../Parallel.hh"
#include "ModelBase.hh"
#include <utility>

//FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
//Model is given a "traits" class that contains stencil, data, force and boundary information

template<class lattice>
using DefaultTraitFlowField = typename DefaultTrait<lattice> :: template SetBoundary<BounceBack>;

template<class lattice, class traits = DefaultTraitFlowField<lattice>>
class FlowField : public CollisionBase<lattice,typename traits::Stencil>, public ModelBase<lattice, traits> { //Inherit from base class to avoid repetition of common
                                                         //calculations
                                                         
    public:

        inline void collide() override; //Collision step

        inline void initialise() override; //Initialisation step

        inline void computeMomenta() override; //Momenta (density, velocity) calculation

        inline const double& getDensity(const int k) const; //Return density at lattice point k

        inline const std::vector<double>& getVelocity() const; //Return vector of velocity

        template<class,class>
        friend class FlowFieldBinary;

    private:

        inline double computeEquilibrium(const double& density, const double* velocity,
                                  const int idx, const int k) const; //Calculate equilibrium in direction idx with a given
                                                        //density and velocity


        inline double computeDensity(const double* distribution, const int k) const; //Calculate density

        inline double computeVelocity(const double* distribution, const double& density,
                                const int xyz, const int k) const; //Calculate velocity

        static constexpr double m_Tau = 1.0; //TEMPORARY relaxation time
        static constexpr double m_InverseTau = 1.0 / m_Tau; //TEMPORARY inverse relaxation time
    
        std::vector<double>& density = Density<>::get<typename traits::Lattice>(); //Reference to vector of densities
        std::vector<double>& velocity = Velocity<>::get<typename traits::Lattice,traits::Lattice::m_NDIM>(); //Reference to vector of velocities

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions
        
};


template<class lattice, class traits>
inline const double& FlowField<lattice, traits>::getDensity(const int k) const {

    return density[k]; //Return reference to density at point k

}

template<class lattice, class traits>
inline const std::vector<double>& FlowField<lattice, traits>::getVelocity() const {

    return velocity; //Return reference to velocity vector

}

template<class lattice, class traits>
inline void FlowField<lattice, traits>::collide() { //Collision step

    #pragma omp for schedule(guided)
    for (int k = lattice::m_HaloSize; k < lattice::m_N - lattice::m_HaloSize; k++) { //loop over k

        if(!ModelBase<lattice, traits>::m_Geometry.isSolid(k)){

            auto forcemethods=ModelBase<lattice,traits>::getForceCalculator(ModelBase<lattice,traits>::mt_Forces,k);

            double* old_distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionOldPointer(k);

            double equilibriums[traits::Stencil::Q];
            double forces[traits::Stencil::Q];

            for (int idx = 0; idx <traits::Stencil::Q; idx++) {
                equilibriums[idx]=computeEquilibrium(density[k], &velocity[k*traits::Lattice::m_NDIM], idx, k);
                if constexpr(std::tuple_size<typename std::remove_reference<decltype(*forcemethods)>::type>::value != 0){
                    forces[idx]=std::apply([idx, k](auto&... forcetype){
                                return (forcetype.template compute<traits>(idx, k) + ...);
                            }, *forcemethods);
                }
            }
            
            this->collisionQ(forces,equilibriums,old_distribution,m_InverseTau,k);

        }
        
    }

    ModelBase<lattice, traits>::m_Data.communicateDistribution();

}

template<class lattice, class traits>
inline void FlowField<lattice, traits>::initialise() { //Initialise model

    ModelBase<lattice, traits>::m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    traits::template CollisionModel<typename traits::Stencil>::template initialise<typename traits::Lattice>(m_InverseTau,m_InverseTau);

    #pragma omp parallel for schedule(guided)
    for (int k = lattice::m_HaloSize; k <lattice::m_N - lattice::m_HaloSize; k++) { //loop over k

        double* distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionOldPointer(k);

        Density<>::initialise<lattice>(1.0,k); //Set density to 1 initially (This will change)
        Velocity<>::initialise<lattice,lattice::m_NDIM>(0.0,k,x);
        Velocity<>::initialise<lattice,lattice::m_NDIM>(0.0,k,y);
        if constexpr (lattice::m_NDIM==3) Velocity<>::initialise<lattice,lattice::m_NDIM>(0.0,k,z);

        for (int idx = 0; idx <traits::Stencil::Q; idx++) {

            double equilibrium = computeEquilibrium(density[k], &velocity[k * traits::Stencil::D], idx, k);

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;        

        }
        
    }
    
}


template<class lattice, class traits>
inline void FlowField<lattice, traits>::computeMomenta() { //Calculate Density<> and Velocity

    #pragma omp for schedule(guided)
    for (int k = lattice::m_HaloSize; k <lattice::m_N - lattice::m_HaloSize; k++) { //Loop over k

        if(!Geometry<lattice>::isSolid(k)){

            double* distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionPointer(k);
            velocity[k * traits::Stencil::D + x] = computeVelocity(distribution, density[k], x, k); //Calculate velocities
            velocity[k * traits::Stencil::D + y] = computeVelocity(distribution,density[k], y, k);
            if constexpr (lattice::m_NDIM == 3) velocity[k * traits::Stencil::D + z] = computeVelocity(distribution ,density[k], z, k);
            density[k] = computeDensity(distribution, k); //Calculate density
            
            
        }

    }

    //ModelBase<lattice, traits>::m_Data.communicate(Density<>::getInstance<lattice>());

}


template<class lattice, class traits>
inline double FlowField<lattice, traits>::computeEquilibrium(const double& density, const double* velocity, const int idx, const int k) const {

    return density * CollisionBase<lattice,typename traits::Stencil>::computeGamma(velocity, idx); //Equilibrium is density
                                                                                        //times gamma in this
                                                                                        //case

}

template<class lattice, class traits>
inline double FlowField<lattice, traits>::computeDensity(const double* distribution, const int k) const { //Density<> calculation
    //Density<> is the sum of distributions plus any source/correction terms

    if constexpr(std::tuple_size<typename traits::Forces>::value != 0) {

        return CollisionBase<lattice,typename traits::Stencil>::computeZerothMoment(distribution) + std::apply([k](auto&... forces) {

                return (forces.template computeDensitySource<traits>(k) + ...);

            }, ModelBase<lattice, traits>::mt_Forces);

    }
    else return CollisionBase<lattice,typename traits::Stencil>::computeZerothMoment(distribution);

}

template<class lattice, class traits>
inline double FlowField<lattice, traits>::computeVelocity(const double* distribution, const double& density,
                                             const int xyz, const int k) const { //Velocity calculation in direction xyz
    //Velocity in direction xyz is sum of distribution times the xyz component of the discrete velocity vector
    //in each direction plus any source/correction terms

    if constexpr(std::tuple_size<typename traits::Forces>::value != 0) {

        return (1./(Density<>::get<lattice>(k)))*CollisionBase<lattice,typename traits::Stencil>::computeFirstMoment(distribution, xyz) + (1./(Density<>::get<lattice>(k)))*std::apply([xyz, k](auto&&... forces) {

                return (forces.template computeVelocitySource<traits>(xyz, k) + ...);
                
            }, ModelBase<lattice, traits>::mt_Forces);

    }
    else return (1./(Density<>::get<lattice>(k)))*CollisionBase<lattice,typename traits::Stencil>::computeFirstMoment(distribution, xyz);

}