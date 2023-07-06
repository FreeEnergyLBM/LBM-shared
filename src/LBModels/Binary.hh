#pragma once
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include "../Parallel.hh"
#include "ModelBase.hh"
#include "../BoundaryModels/Boundaries.hh"
#include "../AddOns/AddOns.hh"
#include "../Forces/Forces.hh"
#include "../GradientStencils/GradientStencils.hh"
#include <utility>
#include <array>
#include <omp.h>


//Binary.hh: Contains the details of the LBM model to solve an equation for phase separation. Each
//Model is given a "T_traits" class that contains stencil, data, force and boundary information

template<class T_lattice>
using DefaultTraitBinary = typename DefaultTrait<T_lattice,2> :: template SetBoundary<BounceBack> 
                                                          :: template SetPreProcessor<ChemicalPotentialCalculatorBinary,CubicWetting> 
                                                          :: template SetPostProcessor<GradientsMultiStencil<OrderParameter<>,CentralXYZ,LaplacianCentral>>;


template<class T_lattice, class T_traits = DefaultTraitBinary<T_lattice>>
class Binary: public CollisionBase<T_lattice, typename T_traits::Stencil>, public ModelBase<T_lattice, T_traits> { //Inherit from base class to avoid repetition of common
                                                      //calculations

    using Stencil = typename T_traits::Stencil;  
    static constexpr int m_NDIM = T_lattice::NDIM; 
                                                      
    public:

        inline void setTau1(double val){m_Tau1=val;}
        inline void setTau2(double val){m_Tau2=val;}

        inline void collide() override; //Collision step

        inline void initialise() override; //Initialisation step

        inline void computeMomenta() override; //Momenta (density, velocity) calculation

    private:

        inline double computeEquilibrium(const double& orderparam, const double* velocity, int idx, int k); //Calculate equilibrium in direction idx with a given density and velocity

        inline double computeModelForce(int xyz, int k); //Calculate forces specific to the model in direction xyz

        inline double computeCollisionQ(const double* forcemethods,const double* equilibriums, int k, const double* old, int idx); //Calculate collision
                                                                                           //at index idx

        static constexpr double m_Tau = 1.0; //TEMPORARY relaxation time
        static constexpr double m_InverseTau = 1.0 / m_Tau; //TEMPORARY inverse relaxation time

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions
        
        double m_Gamma = 1;

        double m_Tau1=1;
        double m_Tau2=1;

        std::vector<double>& orderparameter = OrderParameter<>::get<T_lattice>(); //Reference to vector of order parameters
        std::vector<double>& velocity = Velocity<>::get<T_lattice,T_lattice::NDIM>(); //Reference to vector of velocities
        std::vector<double>& itau = InverseTau<>::get<T_lattice>(); //Reference to vector of velocities

};

template<class T_lattice, class T_traits>
inline void Binary<T_lattice, T_traits>::collide() {

    #pragma omp for schedule(guided)
    for (int k = T_lattice::HaloSize; k < T_lattice::N - T_lattice::HaloSize; k++){ //loop over k
        
        if(!Geometry<T_lattice>::isSolid(k)){
            
            auto forcemethods = this -> getForceCalculator(this -> mt_Forces, k);

            double* old_distribution = this -> m_Distribution.getDistributionOldPointer(k);

            double equilibriumsum = 0;

            double equilibriums[Stencil::Q] = {};
            double forces[Stencil::Q] = {};

            for (int idx = 1; idx <Stencil::Q; idx++) {

                equilibriums[idx] = computeEquilibrium(orderparameter[k], &velocity[k * Stencil::D], idx, k);
                equilibriumsum += equilibriums[idx];

                this -> updateForces(forces[idx], *forcemethods, k, idx);

            }

            equilibriums[0] = orderparameter[k] - equilibriumsum;

            this -> updateForces(forces[0], *forcemethods, k, 0);

            this -> collisionQ(forces, equilibriums, old_distribution, m_InverseTau, k);

        }
        
    }

    this -> m_Data.communicateDistribution();

}

template<class T_lattice, class T_traits>
inline void Binary<T_lattice, T_traits>::initialise() { //Initialise model

    this -> m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    T_traits::template CollisionModel<Stencil>::template initialise<T_lattice>(m_Tau1, m_Tau2);
    
    #pragma omp parallel for schedule(guided) 
    for (int k = T_lattice::HaloSize; k < T_lattice::N - T_lattice::HaloSize; k++) { //loop over k

        double* distribution = this -> m_Distribution.getDistributionPointer(k);
        double* old_distribution = this -> m_Distribution.getDistributionOldPointer(k);

        ChemicalPotential<>::initialise<T_lattice>(0,k);

        OrderParameter<>::initialise<T_lattice>(1.0,k);
        
        InverseTau<>::initialise<T_lattice>(1.0 / (0.5 * (1.0 + OrderParameter<>::get<T_lattice>(k)) * (m_Tau1 - m_Tau2) + m_Tau2), k);

        double equilibriumsum = 0;

        for (int idx = Stencil::Q - 1; idx >= 0; idx--) {

            double equilibrium;

            if (idx> 0) equilibrium = computeEquilibrium(orderparameter[k], &velocity[k * Stencil::D], idx, k);
            else equilibrium = orderparameter[k] - equilibriumsum;

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;        

            equilibriumsum += equilibrium;

        }
        
    }

    this -> m_Data.communicate(SolidLabels<>::getInstance<T_lattice>());

}


template<class T_lattice, class T_traits>
inline void Binary<T_lattice, T_traits>::computeMomenta() { //Calculate order parameter

    #pragma omp for schedule(guided)
    for (int k = T_lattice::HaloSize; k < T_lattice::N - T_lattice::HaloSize; k++) { //Loop over k

        if(!Geometry<T_lattice>::isSolid(k)){

            double* distribution = this -> m_Distribution.getDistributionPointer(k);

            orderparameter[k] = this -> computeDensity(distribution, k);

            itau[k] = 1.0 / (0.5 * (1.0 + orderparameter[k]) * (m_Tau1) - 0.5 * (-1.0 + orderparameter[k]) * m_Tau2);
        
        }

    }

    this -> m_Data.communicate(OrderParameter<>::getInstance<T_lattice>());
}

template<class T_lattice, class T_traits>
inline double Binary<T_lattice, T_traits>::computeEquilibrium(const double& orderparam, const double* velocity, int idx, int k) {

    return Stencil::Weights[idx] * (ChemicalPotential<>::get<T_lattice>(k) * m_Gamma / Stencil::Cs2 + orderparam * CollisionBase<T_lattice, Stencil>::computeVelocityFactor(velocity, idx));

}

//FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
//Model is given a "T_traits" class that contains stencil, data, force and boundary information

template<class T_lattice>
using DefaultTraitFlowFieldBinary = typename DefaultTrait<T_lattice,2> :: template SetBoundary<BounceBack>
                                                                     :: template SetForce<ChemicalForce<Guo,Gradient>>;

template<class T_lattice, class T_traits = DefaultTraitFlowFieldBinary<T_lattice>>
class FlowFieldBinary : public FlowField<T_lattice, T_traits>{ //Inherit from base class to avoid repetition of common
                                                         //calculations

    using Stencil = typename T_traits::Stencil;  
    static constexpr int m_NDIM = T_lattice::NDIM; 
    
    public:

        inline void setTau1(double val){m_Tau1 = val;}
        inline void setTau2(double val){m_Tau2 = val;}

        inline virtual void collide() override; //Collision step

        inline virtual void initialise() override; //Initialisation step

    private:

        inline double computeEquilibrium(const double& density, const double* velocity, const double& order_parameter, const double& chemical_potential, int idx, int k); //Calculate equilibrium in direction idx with a given//density and velocity
                                                                     //at index idx

        double m_Tau1 = 1;
        double m_Tau2 = 1;

        enum{x = 0, y = 1, z = 2};

};

template<class T_lattice, class T_traits>
inline double FlowFieldBinary<T_lattice, T_traits>::computeEquilibrium(const double& density, const double* velocity, const double& order_parameter, const double& chemical_potential, int idx, int k) {

    return density * CollisionBase<T_lattice, Stencil>::computeGamma(velocity, idx) + Stencil::Weights[idx] * order_parameter * chemical_potential / Stencil::Cs2; //Equilibrium is density times gamma in this case

}

template<class T_lattice, class T_traits>
inline void FlowFieldBinary<T_lattice, T_traits>::collide() { //Collision step
       
    #pragma omp for schedule(guided)
    for (int k = T_lattice::HaloSize; k <T_lattice::N - T_lattice::HaloSize; k++) { //loop over k

        if(!Geometry<T_lattice>::isSolid(k)){

            auto forcemethods=this -> getForceCalculator(this -> mt_Forces,k);

            //int QQ = Stencil::Q;
            double* old_distribution = this -> m_Distribution.getDistributionOldPointer(k);
            double equilibriumsum = 0;
            
            double equilibriums[Stencil::Q];
            double forces[Stencil::Q];

            for (int idx = 1; idx < Stencil::Q; idx++) {

                equilibriums[idx] = computeEquilibrium(this -> density[k], &(this -> velocity[k * Stencil::D]), OrderParameter<>::get<T_lattice>(k), ChemicalPotential<>::get<T_lattice>(k), idx, k);
                equilibriumsum += equilibriums[idx];

                this -> updateForces(forces[idx], *forcemethods, k, idx);

            }

            this -> updateForces(forces[0], *forcemethods, k, 0);

            equilibriums[0]=this -> density[k]-equilibriumsum;

            this -> collisionQ(forces,equilibriums,old_distribution,InverseTau<>::get<T_lattice>(k),k);

        }
        
    }

    FlowField<T_lattice, T_traits>::m_Data.communicateDistribution();

}

template<class T_lattice, class T_traits>
inline void FlowFieldBinary<T_lattice, T_traits>::initialise() { //Initialise model

    this -> m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    T_traits::template CollisionModel<Stencil>::template initialise<T_lattice>(m_Tau1,m_Tau2);

    #pragma omp parallel for schedule(guided)
    for (int k = T_lattice::HaloSize; k <T_lattice::N - T_lattice::HaloSize; k++) { //loop over k

        double* distribution = this -> m_Distribution.getDistributionPointer(k);
        double* old_distribution = this -> m_Distribution.getDistributionOldPointer(k);

        Density<>::initialise<T_lattice>(1.0, k); //Set density to 1 initially (This will change)
        Velocity<>::initialise<T_lattice,T_lattice::NDIM>(0.0, k, x);
        Velocity<>::initialise<T_lattice,T_lattice::NDIM>(0.0, k, y);
        if constexpr (T_lattice::NDIM == 3) Velocity<>::initialise<T_lattice, T_lattice::NDIM>(0.0, k, z);

        double equilibriumsum = 0;

        for (int idx = Stencil::Q-1; idx>= 0; idx--) {

            double equilibrium;

            if (idx>0) equilibrium = computeEquilibrium(this -> density[k], &(this -> velocity[k * Stencil::D]), OrderParameter<>::get<T_lattice>(k), ChemicalPotential<>::get<T_lattice>(k), idx, k);

            else equilibrium = FlowField<T_lattice, T_traits>::density[k] - equilibriumsum;

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;
            
            equilibriumsum += equilibrium;
        }
        
    }
    
}
