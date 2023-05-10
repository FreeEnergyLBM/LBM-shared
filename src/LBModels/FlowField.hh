#pragma once
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include "../BoundaryModels/Boundaries.hh"
#include "../Forces/Forces.hh"
#include "../Parallel.hh"
#include "ModelBase.hh"
#include <utility>

//FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
//Model is given a "traits" class that contains stencil, data, force and boundary information

template<class properties>
struct DefaultTraitFlowField{
    
    using Stencil = std::conditional_t<std::remove_reference<properties>::type::m_NDIM == 2, D2Q9, D3Q19>; //Here, D refers to the number of cartesian dimensions

    using Boundaries = std::tuple<BounceBack>;

    using Forces = std::tuple<>;

    using Properties = properties;

};



template<class traits = DefaultTraitFlowField<decltype(GETPROPERTIES())>>
class FlowField : public CollisionBase<typename traits::Stencil>, public ModelBase<traits> { //Inherit from base class to avoid repetition of common
                                                         //calculations
                                                         
    public:

        void collide() override; //Collision step

        void initialise() override; //Initialisation step

        void computeMomenta(); //Momenta (density, velocity) calculation

        const double& getDensity(const int k) const; //Return density at lattice point k

        const std::vector<double>& getVelocity() const; //Return vector of velocity

        template<class>
        friend class FlowFieldBinary;

    private:

       double computeEquilibrium(const double& density, const double* velocity,
                                  const int idx, const int k) const; //Calculate equilibrium in direction idx with a given
                                                        //density and velocity

        double computeModelForce(int xyz, int k) const; //Calculate forces specific to the model in direction xyz

        double computeCollisionQ(int k, const double& old, const double& density,
                                 const double* velocity, const int idx) const; //Calculate collision
                                                                                           //at index idx

        double computeDensity(const double* distribution, const int k) const; //Calculate density

        double computeVelocity(const double* distribution, const double& density,
                                const int xyz, const int k) const; //Calculate velocity

        static constexpr double m_Tau = 1.0; //TEMPORARY relaxation time
        static constexpr double m_InverseTau = 1.0 / m_Tau; //TEMPORARY inverse relaxation time

        Density m_Density; //Density<>
        Velocity m_Velocity; //Velocity
    
        vector<double>& density = m_Density.getParameter(); //Reference to vector of densities
        vector<double>& velocity = m_Velocity.getParameter(); //Reference to vector of velocities

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions
        
};


template<class traits>
const double& FlowField<traits>::getDensity(const int k) const {

    return density[k]; //Return reference to density at point k

}

template<class traits>
const std::vector<double>& FlowField<traits>::getVelocity() const {

    return velocity; //Return reference to velocity vector

}

template<class traits>
void FlowField<traits>::collide() { //Collision step

    #pragma omp parallel for schedule(guided)
    for (int k = GETPROPERTIES().m_HaloSize; k <GETPROPERTIES().m_N - GETPROPERTIES().m_HaloSize; k++) { //loop over k

        double* old_distribution = ModelBase<traits>::m_Distribution.getDistributionOldPointer(k);

        for (int idx = 0; idx <traits::Stencil::Q; idx++) { //loop over discrete velocity directions
            //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"
            double collision = computeCollisionQ(k, old_distribution[idx], density[k], &velocity[k * traits::Stencil::D], idx);

            ModelBase<traits>::m_Distribution.getDistributionPointer(ModelBase<traits>::m_Distribution.streamIndex(k, idx))[idx] = collision;

        }
        
    }

    #ifdef MPIPARALLEL
    ModelBase<traits>::m_Data.communicateDistribution();
    #endif

}

template<class traits>
void FlowField<traits>::initialise() { //Initialise model

    ModelBase<traits>::m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    
    #pragma omp parallel for schedule(guided)
    for (int k = GETPROPERTIES().m_HaloSize; k <GETPROPERTIES().m_N - GETPROPERTIES().m_HaloSize; k++) { //loop over k

        double* distribution = ModelBase<traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution = ModelBase<traits>::m_Distribution.getDistributionOldPointer(k);

        density[k] = 1.0; //Set density to 1 initially (This will change)
        velocity[k * traits::Stencil::D + x] = 0.0; //0 initial velocity
        velocity[k * traits::Stencil::D + y] = 0.0;
        velocity[k * traits::Stencil::D + z] = 0.0;

        for (int idx = 0; idx <traits::Stencil::Q; idx++) {

            double equilibrium = computeEquilibrium(density[k], &velocity[k * traits::Stencil::D], idx, k);

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;        

        }
        
    }
    
}


template<class traits>
void FlowField<traits>::computeMomenta() { //Calculate Density<> and Velocity

    #pragma omp parallel for schedule(guided)
    for (int k = GETPROPERTIES().m_HaloSize; k <GETPROPERTIES().m_N - GETPROPERTIES().m_HaloSize; k++) { //Loop over k

        double* distribution = ModelBase<traits>::m_Distribution.getDistributionPointer(k);

        density[k] = computeDensity(distribution, k); //Calculate density
        velocity[k * traits::Stencil::D + x] = computeVelocity(distribution, density[k], x, k); //Calculate velocities
        velocity[k * traits::Stencil::D + y]=computeVelocity(distribution,density[k], y, k);
        if constexpr (traits::Stencil::D == 3) velocity[k * traits::Stencil::D + z] = computeVelocity(distribution ,density[k], z, k);

    }

    #ifdef MPIPARALLEL
    ModelBase<traits>::m_Data.communicate(m_Density);
    #endif

}

template<class traits>
double FlowField<traits>::computeCollisionQ(const int k, const double& old, const double& density,
                                               const double* velocity, const int idx) const {
                                            //Calculate collision step at a given velocity index at point k

    double forcexyz[traits::Stencil::D]; //Temporary array storing force in each cartesian direction

    //Force is the sum of model forces and given forces
    for(int xyz = 0; xyz <traits::Stencil::D; xyz++) forcexyz[xyz] = computeModelForce(xyz, k) + ModelBase<traits>::computeForces(xyz, k);
    
    //Sum of collision + force contributions
    return CollisionBase<typename traits::Stencil>::collideSRT(old, computeEquilibrium(density, velocity, idx, k), m_InverseTau)
              + CollisionBase<typename traits::Stencil>::forceGuoSRT(forcexyz, velocity, m_InverseTau, idx);

}


template<class traits>
double FlowField<traits>::computeEquilibrium(const double& density, const double* velocity, const int idx, const int k) const {

    return density * CollisionBase<typename traits::Stencil>::computeGamma(velocity, idx); //Equilibrium is density
                                                                                        //times gamma in this
                                                                                        //case

}

template<class traits>
double FlowField<traits>::computeModelForce(int k, int xyz) const {

    return 0.0; //No model force in this case

}

template<class traits>
double FlowField<traits>::computeDensity(const double* distribution, const int k) const { //Density<> calculation
    //Density<> is the sum of distributions plus any source/correction terms

    if constexpr(std::tuple_size<typename traits::Forces>::value != 0) {

        return CollisionBase<typename traits::Stencil>::computeZerothMoment(distribution) + std::apply([k](auto&... forces) {

                return (forces.computeDensitySource(k) + ...);

            }, ModelBase<traits>::mt_Forces);

    }
    else return CollisionBase<typename traits::Stencil>::computeZerothMoment(distribution);

}

template<class traits>
double FlowField<traits>::computeVelocity(const double* distribution, const double& density,
                                             const int xyz, const int k) const { //Velocity calculation in direction xyz
    //Velocity in direction xyz is sum of distribution times the xyz component of the discrete velocity vector
    //in each direction plus any source/correction terms

    if constexpr(std::tuple_size<typename traits::Forces>::value != 0) {

        return CollisionBase<typename traits::Stencil>::computeFirstMoment(distribution, xyz) + std::apply([xyz, k](auto&&... forces) {

                return (forces.computeVelocitySource(xyz, k) + ...);
                
            }, ModelBase<traits>::mt_Forces);

    }
    else return CollisionBase<typename traits::Stencil>::computeFirstMoment(distribution, xyz);

}