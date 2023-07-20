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

            //auto forcemethods = this -> getForceCalculator(this -> mt_Forces, k);

            double* old_distribution = this -> m_Distribution.getDistributionOldPointer(k);

            double equilibriums[Stencil::Q];
            //double forces[Stencil::Q]; 

            for (int idx = 0; idx < Stencil::Q; idx++) {

                equilibriums[idx] = computeEquilibrium(density[k], &velocity[k * m_NDIM], idx, k);

                //this -> updateForces(forces[idx], *forcemethods, k, idx);
            }
            
            this -> collisionQ(equilibriums, old_distribution, m_InverseTau,k); // CHANGE NEEDED If no forces, don't require them to be passed

        }
        
    }

    this -> m_Data.communicateDistribution();

}

template<class T_lattice, class T_traits>
inline void FlowField<T_lattice, T_traits>::initialise() { //Initialise model

    this -> m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    T_traits::template CollisionModel<Stencil>::template initialise<T_lattice>(this -> mt_Forces,m_Tau,m_Tau);

    #pragma omp parallel for schedule(guided)
    for (int k = T_lattice::HaloSize; k < T_lattice::N - T_lattice::HaloSize; k++) { //loop over k

        double* distribution = this -> m_Distribution.getDistributionPointer(k);
        double* old_distribution = this -> m_Distribution.getDistributionOldPointer(k);

        Density<>::initialise<T_lattice>(1.0 ,k); //Set density to 1 initially (This will change)
        Velocity<>::initialise<T_lattice, m_NDIM>(0.0, k, x);
        if constexpr (m_NDIM >= 2) Velocity<>::initialise<T_lattice, m_NDIM>(0.0, k, y);
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
            
            if constexpr (m_NDIM >= 2) velocity[k * Stencil::D + y] = this -> computeVelocity(distribution,density[k], y, k);
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

template<class method>
class PressureForce : public ChemicalForce<method> {
    public:
        template<class traits>
        inline double computeXYZ(const int xyz, const int k) {
            
            return traits::Stencil::Cs2*GradientDensity<>::get<typename traits::Lattice, traits::Lattice::NDIM>(k,xyz)+ChemicalForce<method>::template computeXYZ<traits>(xyz,k);
        
        }
        template<class traits>
        inline double computeQ(const int idx, const int k) {
            double sum=0;
            for(int xyz=0;xyz<traits::Lattice::NDIM;xyz++){
                sum+=(traits::Stencil::Cs2*GradientDensity<>::get<typename traits::Lattice, traits::Lattice::NDIM>(k,xyz)+ChemicalForce<method>::template computeXYZ<traits>(xyz,k))*traits::Stencil::Ci_xyz(xyz)[idx];
            }
            return sum;
        
        }
        template<class traits>
        inline double computeDensitySource(int k) { //SHOULD BE CENTRAL GRADIENTS
            double source=0;
            for(int xyz = 0; xyz<traits::Lattice::NDIM; xyz++) source+=traits::Lattice::DT*0.5*traits::Stencil::Cs2*GradientDensity<>::get<typename traits::Lattice, traits::Lattice::NDIM>(k,xyz)*Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,xyz);
            return source;
        }

    private:

};

template<class lattice>
using DefaultTraitFlowFieldPressure = typename DefaultTrait<lattice> ::template AddPreProcessor<Gradients<Density<>,CentralXYZ>> ::template AddForce<PressureForce<He>>;

template<class lattice, class traits = DefaultTraitFlowFieldPressure<lattice>>
class FlowFieldPressure : public CollisionBase<lattice,typename traits::Stencil>, public ModelBase<lattice, traits> { //Inherit from base class to avoid repetition of common
                                                         //calculations

    using Stencil = typename traits::Stencil;  
    static constexpr int m_NDIM = lattice::NDIM; 

    public:

        inline void setTau1(double val){m_Tau1 = val;}
        inline void setTau2(double val){m_Tau2 = val;}

        inline void collide() override; //Collision step

        inline void initialise() override; //Initialisation step

        inline void computeMomenta() override; //Momenta (density, velocity) calculation

        inline double computeEquilibrium(const double& density, const double& pressure, const double* velocity,
                                  const int idx, const int k) const; //Calculate equilibrium in direction idx with a given
                                                        //density and velocitys

        std::vector<double>& pressure = Pressure<>::get<lattice>(); //Reference to vector of densities
        std::vector<double>& density = Density<>::get<lattice>(); //Reference to vector of densities
        std::vector<double>& velocity = Velocity<>::get<lattice,lattice::NDIM>(); //Reference to vector of velocities

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions

    private:
        double m_Tau1 = 1;
        double m_Tau2 = 1;
        
};

template<class lattice, class traits>
inline void FlowFieldPressure<lattice, traits>::collide() { //Collision step

    #pragma omp for schedule(guided)
    for (int k = lattice::HaloSize; k < lattice::N - lattice::HaloSize; k++) { //loop over k

        if(!ModelBase<lattice, traits>::m_Geometry.isSolid(k)){

            double* old_distribution = this -> m_Distribution.getDistributionOldPointer(k);

            double equilibriums[Stencil::Q];

            for (int idx = 0; idx < Stencil::Q; idx++) {

                equilibriums[idx] = computeEquilibrium(density[k], pressure[k], &velocity[k * m_NDIM], idx, k);

            }
            
            this -> collisionQ(equilibriums, old_distribution, InverseTau<>::get<lattice>(k),k); // CHANGE NEEDED If no forces, don't require them to be passed

        }
        
    }

    ModelBase<lattice, traits>::m_Data.communicateDistribution();

}

template<class lattice, class traits>
inline void FlowFieldPressure<lattice, traits>::initialise() { //Initialise model

    ModelBase<lattice, traits>::m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    traits::template CollisionModel<Stencil>::template initialise<lattice>(this -> mt_Forces,m_Tau1,m_Tau2);
    
    #pragma omp parallel for schedule(guided)
    for (int k = 0; k <lattice::N; k++) { //loop over k

        double* distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionOldPointer(k);
        InverseTau<>::initialise<lattice>(1.0,k);
        Pressure<>::initialise<lattice>(1.0,k); //Set density to 1 initially (This will change)
        Density<>::initialise<lattice>(1.0,k);
        Velocity<>::initialise<lattice,lattice::NDIM>(0.0,k,x);
        Velocity<>::initialise<lattice,lattice::NDIM>(0.0,k,y);
        if constexpr (lattice::NDIM==3) Velocity<>::initialise<lattice,lattice::NDIM>(0.0,k,z);

        for (int idx = 0; idx <traits::Stencil::Q; idx++) {

            double equilibrium = computeEquilibrium(density[k], pressure[k], &velocity[k * traits::Stencil::D], idx, k);

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;        

        }
        
    }

    ModelBase<lattice, traits>::m_Data.communicate(Pressure<>::get<lattice>());
    ModelBase<lattice, traits>::m_Data.communicate(Velocity<>::get<lattice,lattice::NDIM>());
    
}


template<class lattice, class traits>
inline void FlowFieldPressure<lattice, traits>::computeMomenta() { //Calculate Density<> and Velocity

    #pragma omp for schedule(guided)
    for (int k = lattice::HaloSize; k <lattice::N - lattice::HaloSize; k++) { //Loop over k

        if(!ModelBase<lattice, traits>::m_Geometry.isSolid(k)){

            double* distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionPointer(k);

            pressure[k] = this -> computeDensity(distribution, k); //Calculate density
            
            velocity[k * traits::Stencil::D + x] = 1./(traits::Stencil::Cs2)*this->computeVelocity(distribution, density[k], x, k); //Calculate velocities
            velocity[k * traits::Stencil::D + y] = 1./(traits::Stencil::Cs2)*this->computeVelocity(distribution,density[k], y, k);
            if constexpr (lattice::NDIM == 3) velocity[k * traits::Stencil::D + z] = 1./(traits::Stencil::Cs2)*this->computeVelocity(distribution ,density[k], z, k);

        }

    }

    ModelBase<lattice, traits>::m_Data.communicate(Pressure<>::get<lattice>());
    ModelBase<lattice, traits>::m_Data.communicate(Velocity<>::get<lattice,lattice::NDIM>());

}

template<class lattice, class traits>
inline double FlowFieldPressure<lattice, traits>::computeEquilibrium(const double& density, const double& pressure, const double* velocity, const int idx, const int k) const {
    
    return traits::Stencil::Weights[idx]*(pressure + density * traits::Stencil::Cs2 * CollisionBase<lattice,typename traits::Stencil>::computeVelocityFactor(velocity, idx)); //Equilibrium is density //Equilibrium is density
                                                                                        //times gamma in this
                                                                                        //case

}