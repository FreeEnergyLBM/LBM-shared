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

template<class T_lattice>
using DefaultTraitFlowField = typename DefaultTrait<T_lattice> :: template SetBoundary<BounceBack>;

template<class T_lattice, class T_traits = DefaultTraitFlowField<T_lattice>>
class FlowField : public CollisionBase<T_lattice,typename T_traits::Stencil>, public ModelBase<T_lattice, T_traits> { //Inherit from base class to avoid repetition of common
                                                         //calculations
                                                         
    using Stencil = typename T_traits::Stencil;  
    static constexpr int m_NDIM = T_lattice::NDIM; 

    public:

        inline void collide() override; //Collision step

        inline void initialise() override; //Initialisation step

        inline void computeMomenta() override; //Momenta (density, velocity) calculation

        template<class,class>
        friend class FlowFieldBinary;

    private:

        inline double computeEquilibrium(const double& density, const double* velocity,
                                  int idx, int k); //Calculate equilibrium in direction idx with a given
                                                        //density and velocity

        static constexpr double m_Tau = 1.0; //TEMPORARY relaxation time
        static constexpr double m_InverseTau = 1.0 / m_Tau; //TEMPORARY inverse relaxation time
    
        std::vector<double>& density = Density<>::get<T_lattice>(); //Reference to vector of densities
        std::vector<double>& velocity = Velocity<>::get<T_lattice,m_NDIM>(); //Reference to vector of velocities

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions
        
};

template<class T_lattice, class T_traits>
inline void FlowField<T_lattice, T_traits>::collide() { //Collision step

    #pragma omp for schedule(guided)
    for (int k = T_lattice::HaloSize; k < T_lattice::N - T_lattice::HaloSize; k++) { //loop over k

        if(!Geometry<T_lattice>::isSolid(k)){

            auto forcemethods = this -> getForceCalculator(this -> mt_Forces, k);

            double* old_distribution = this -> m_Distribution.getDistributionOldPointer(k);

            double equilibriums[Stencil::Q];
            double forces[Stencil::Q]; 

            for (int idx = 0; idx < Stencil::Q; idx++) {

                equilibriums[idx] = computeEquilibrium(density[k], &velocity[k * m_NDIM], idx, k);

                this -> updateForces(forces[idx], *forcemethods, k, idx);
            }
            
            this -> collisionQ(forces, equilibriums, old_distribution, m_InverseTau,k); // CHANGE NEEDED If no forces, don't require them to be passed

        }
        
    }

    this -> m_Data.communicateDistribution();

}

template<class T_lattice, class T_traits>
inline void FlowField<T_lattice, T_traits>::initialise() { //Initialise model

    this -> m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    T_traits::template CollisionModel<Stencil>::template initialise<T_lattice>(m_InverseTau, m_InverseTau);

    #pragma omp parallel for schedule(guided)
    for (int k = T_lattice::HaloSize; k < T_lattice::N - T_lattice::HaloSize; k++) { //loop over k

        double* distribution = this -> m_Distribution.getDistributionPointer(k);
        double* old_distribution = this -> m_Distribution.getDistributionOldPointer(k);

        Density<>::initialise<T_lattice>(1.0 ,k); //Set density to 1 initially (This will change)
        Velocity<>::initialise<T_lattice, m_NDIM>(0.0, k, x);
        Velocity<>::initialise<T_lattice, m_NDIM>(0.0, k, y);
        if constexpr (m_NDIM == 3) Velocity<>::initialise<T_lattice, m_NDIM>(0.0, k, z);

        for (int idx = 0; idx <Stencil::Q; idx++) {

            double equilibrium = computeEquilibrium(density[k], &velocity[k * Stencil::D], idx, k);

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;        

        }
        
    }
    
}


template<class T_lattice, class T_traits>
inline void FlowField<T_lattice, T_traits>::computeMomenta() { //Calculate Density<> and Velocity

    #pragma omp for schedule(guided)
    for (int k = T_lattice::HaloSize; k <T_lattice::N - T_lattice::HaloSize; k++) { //Loop over k

        if(!Geometry<T_lattice>::isSolid(k)){

            double* distribution = this -> m_Distribution.getDistributionPointer(k);
            velocity[k * Stencil::D + x] = this -> computeVelocity(distribution, density[k], x, k); //Calculate velocities
            velocity[k * Stencil::D + y] = this -> computeVelocity(distribution,density[k], y, k);
            if constexpr (m_NDIM == 3) velocity[k * Stencil::D + z] = this -> computeVelocity(distribution ,density[k], z, k);
            density[k] = this -> computeDensity(distribution, k); //Calculate density
            
            
        }

    }

    //this -> m_Data.communicate(Density<>::getInstance<T_lattice>());

}


template<class T_lattice, class T_traits>
inline double FlowField<T_lattice, T_traits>::computeEquilibrium(const double& density, const double* velocity, int idx, int k) {

    return density * CollisionBase<T_lattice,Stencil>::computeGamma(velocity, idx); //Equilibrium is density
                                                                                        //times gamma in this
                                                                                        //case

}