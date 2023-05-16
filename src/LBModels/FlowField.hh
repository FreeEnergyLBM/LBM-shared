#pragma once
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include "../BoundaryModels/Boundaries.hh"
#include "../AddOns/AddOns.hh"
#include "../Parallel.hh"
#include "ModelBase.hh"
#include <utility>

//FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
//Model is given a "traits" class that contains stencil, data, force and boundary information

template<class lattice>
struct DefaultTraitFlowField : BaseTrait<DefaultTraitFlowField<lattice>> {
    
    using Stencil = std::conditional_t<lattice::m_NDIM == 2, D2Q9, D3Q19>; //Here, D refers to the number of cartesian dimensions

    using Boundaries = std::tuple<BounceBack<lattice>>;

    using AddOns = std::tuple<>;

};



template<class lattice, class traits = DefaultTraitFlowField<lattice>>
class FlowField : public CollisionBase<lattice,typename traits::Stencil>, public ModelBase<lattice, traits> { //Inherit from base class to avoid repetition of common
                                                         //calculations
                                                         
    public:

        inline void collide() override; //Collision step

        inline void initialise() override; //Initialisation step

        inline void computeMomenta(); //Momenta (density, velocity) calculation

        inline const double& getDensity(const int k) const; //Return density at lattice point k

        inline const std::vector<double>& getVelocity() const; //Return vector of velocity

        template<class,class>
        friend class FlowFieldBinary;

    private:

        inline double computeEquilibrium(const double& density, const double* velocity,
                                  const int idx, const int k) const; //Calculate equilibrium in direction idx with a given
                                                        //density and velocity

        inline double computeModelForce(int xyz, int k) const; //Calculate forces specific to the model in direction xyz

        inline double computeCollisionQ(int k, const double& old, const double& density,
                                 const double* velocity, const int idx) const; //Calculate collision
                                                                                           //at index idx

        inline double computeDensity(const double* distribution, const int k) const; //Calculate density

        inline double computeVelocity(const double* distribution, const double& density,
                                const int xyz, const int k) const; //Calculate velocity

        static constexpr double m_Tau = 1.0; //TEMPORARY relaxation time
        static constexpr double m_InverseTau = 1.0 / m_Tau; //TEMPORARY inverse relaxation time

        Density<lattice> m_Density; //Density<>
        Velocity<lattice> m_Velocity; //Velocity
    
        std::vector<double>& density = m_Density.getParameter(); //Reference to vector of densities
        std::vector<double>& velocity = m_Velocity.getParameter(); //Reference to vector of velocities

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

            double* old_distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionOldPointer(k);
            
            for (int idx = 0; idx <traits::Stencil::Q; idx++) { //loop over discrete velocity directions
                //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
                //"computeCollisionQ"
                double collision = computeCollisionQ(k, old_distribution[idx], density[k], &velocity[k * traits::Stencil::D], idx);

                ModelBase<lattice, traits>::m_Distribution.getDistributionPointer(ModelBase<lattice, traits>::m_Distribution.streamIndex(k, idx))[idx] = collision;

            }

        }
        
    }

    ModelBase<lattice, traits>::m_Data.communicateDistribution();

}

template<class lattice, class traits>
inline void FlowField<lattice, traits>::initialise() { //Initialise model

    ModelBase<lattice, traits>::m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    
    #pragma omp parallel for schedule(guided)
    for (int k = lattice::m_HaloSize; k <lattice::m_N - lattice::m_HaloSize; k++) { //loop over k

        double* distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionOldPointer(k);

        m_Density.initialise(1.0,k); //Set density to 1 initially (This will change)
        m_Velocity.initialise(0.0,k,x);
        m_Velocity.initialise(0.0,k,y);
        if constexpr (lattice::m_NDIM==3) m_Velocity.initialise(0.0,k,z);

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

        if(!ModelBase<lattice, traits>::m_Geometry.isSolid(k)){

            double* distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionPointer(k);

            density[k] = computeDensity(distribution, k); //Calculate density
            velocity[k * traits::Stencil::D + x] = computeVelocity(distribution, density[k], x, k); //Calculate velocities
            velocity[k * traits::Stencil::D + y] = computeVelocity(distribution,density[k], y, k);
            if constexpr (lattice::m_NDIM == 3) velocity[k * traits::Stencil::D + z] = computeVelocity(distribution ,density[k], z, k);

        }

    }

    ModelBase<lattice, traits>::m_Data.communicate(m_Density);

}

template<class lattice, class traits>
inline double FlowField<lattice, traits>::computeCollisionQ(const int k, const double& old, const double& density,
                                               const double* velocity, const int idx) const {
                                            //Calculate collision step at a given velocity index at point k

    double forcexyz[traits::Stencil::D]; //Temporary array storing force in each cartesian direction

    //Force is the sum of model forces and given forces
    for(int xyz = 0; xyz <traits::Stencil::D; xyz++) forcexyz[xyz] = computeModelForce(xyz, k) + ModelBase<lattice, traits>::computeForces(xyz, k);
    
    //Sum of collision + force contributions
    return CollisionBase<lattice,typename traits::Stencil>::collideSRT(old, computeEquilibrium(density, velocity, idx, k), m_InverseTau)
              + CollisionBase<lattice,typename traits::Stencil>::forceGuoSRT(forcexyz, velocity, m_InverseTau, idx);

}


template<class lattice, class traits>
inline double FlowField<lattice, traits>::computeEquilibrium(const double& density, const double* velocity, const int idx, const int k) const {

    return density * CollisionBase<lattice,typename traits::Stencil>::computeGamma(velocity, idx); //Equilibrium is density
                                                                                        //times gamma in this
                                                                                        //case

}

template<class lattice, class traits>
inline double FlowField<lattice, traits>::computeModelForce(int k, int xyz) const {

    return 0.0; //No model force in this case

}

template<class lattice, class traits>
inline double FlowField<lattice, traits>::computeDensity(const double* distribution, const int k) const { //Density<> calculation
    //Density<> is the sum of distributions plus any source/correction terms

    if constexpr(std::tuple_size<typename traits::AddOns>::value != 0) {

        return CollisionBase<lattice,typename traits::Stencil>::computeZerothMoment(distribution) + std::apply([k](auto&... addons) {

                return (addons.computeDensitySource(k) + ...);

            }, ModelBase<lattice, traits>::mt_AddOns);

    }
    else return CollisionBase<lattice,typename traits::Stencil>::computeZerothMoment(distribution);

}

template<class lattice, class traits>
inline double FlowField<lattice, traits>::computeVelocity(const double* distribution, const double& density,
                                             const int xyz, const int k) const { //Velocity calculation in direction xyz
    //Velocity in direction xyz is sum of distribution times the xyz component of the discrete velocity vector
    //in each direction plus any source/correction terms

    if constexpr(std::tuple_size<typename traits::AddOns>::value != 0) {

        return CollisionBase<lattice,typename traits::Stencil>::computeFirstMoment(distribution, xyz) + std::apply([xyz, k](auto&&... addons) {

                return (addons.computeVelocitySource(xyz, k) + ...);
                
            }, ModelBase<lattice, traits>::mt_AddOns);

    }
    else return CollisionBase<lattice,typename traits::Stencil>::computeFirstMoment(distribution, xyz);

}
