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
//Model is given a "TTraits" class that contains stencil, data, force and boundary information

template<class TLattice>
using DefaultTraitBinary = typename DefaultTrait<TLattice,2> :: template SetBoundary<BounceBack> 
                                                          :: template SetPreProcessor<ChemicalPotentialCalculatorBinary,CubicWetting> 
                                                          :: template SetPostProcessor<GradientsMultiStencil<OrderParameter<>,CentralXYZ,LaplacianCentral>>;


template<class TLattice, class TTraits = DefaultTraitBinary<TLattice>>
class Binary: public CollisionBase<TLattice, typename TTraits::Stencil>, public ModelBase<TLattice, TTraits> { //Inherit from base class to avoid repetition of common
                                                      //calculations

    using Stencil = typename TTraits::Stencil;  
    static constexpr int mNDIM = TLattice::NDIM; 
                                                      
    public:

        inline void setTau1(double val){mTau1=val;}
        inline void setTau2(double val){mTau2=val;}

        inline void collide() override; //Collision step

        inline void initialise() override; //Initialisation step

        inline void computeMomenta() override; //Momenta (density, velocity) calculation

    private:

        inline double computeEquilibrium(const double& orderparam, const double* velocity, int idx, int k); //Calculate equilibrium in direction idx with a given density and velocity

        inline double computeModelForce(int xyz, int k); //Calculate forces specific to the model in direction xyz

        static constexpr double mTau = 1.0; //TEMPORARY relaxation time
        static constexpr double mInverseTau = 1.0 / mTau; //TEMPORARY inverse relaxation time

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions
        
        double mGamma = 1;

        double mTau1=1;
        double mTau2=1;

        std::vector<double>& orderparameter = OrderParameter<>::get<TLattice>(); //Reference to vector of order parameters
        std::vector<double>& velocity = Velocity<>::get<TLattice,TLattice::NDIM>(); //Reference to vector of velocities
        std::vector<double>& itau = InverseTau<>::get<TLattice>(); //Reference to vector of velocities

};

template<class TLattice, class TTraits>
inline void Binary<TLattice, TTraits>::collide() {

    #pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++){ //loop over k
        
        if(!Geometry<TLattice>::isSolid(k)){

            double* old_distribution = this -> mDistribution.getDistributionOldPointer(k);

            double sum = 0;

            double equilibriums[Stencil::Q] = {};

            for (int idx = 1; idx < Stencil::Q; idx++) {

                equilibriums[idx] = computeEquilibrium(orderparameter[k], &velocity[k * Stencil::D], idx, k);
                sum += equilibriums[idx];

            }

            equilibriums[0] = orderparameter[k] - sum;

            this -> collisionQ(equilibriums, old_distribution, mInverseTau, k);

        }
        
    }

    this -> mData.communicateDistribution();

}

template<class TLattice, class TTraits>
inline void Binary<TLattice, TTraits>::initialise() { //Initialise model

    this -> mData.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this -> mt_Forces,mTau1,mTau2);
    
    #pragma omp parallel for schedule(guided) 
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) { //loop over k

        double* distribution = this -> mDistribution.getDistributionPointer(k);
        double* old_distribution = this -> mDistribution.getDistributionOldPointer(k);

        ChemicalPotential<>::initialise<TLattice>(0,k);

        OrderParameter<>::initialise<TLattice>(1.0,k);
        
        InverseTau<>::initialise<TLattice>(1.0 / (0.5 * (1.0 + OrderParameter<>::get<TLattice>(k)) * (mTau1 - mTau2) + mTau2), k);

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

    this -> mData.communicate(SolidLabels<>::getInstance<TLattice>());

}


template<class TLattice, class TTraits>
inline void Binary<TLattice, TTraits>::computeMomenta() { //Calculate order parameter

    #pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) { //Loop over k

        if(!Geometry<TLattice>::isSolid(k)){

            double* distribution = this -> mDistribution.getDistributionPointer(k);

            orderparameter[k] = this -> computeDensity(distribution, k);

            itau[k] = 1.0 / (0.5 * (1.0 + orderparameter[k]) * (mTau1) - 0.5 * (-1.0 + orderparameter[k]) * mTau2);
        
        }

    }

    this -> mData.communicate(OrderParameter<>::getInstance<TLattice>());
}

template<class TLattice, class TTraits>
inline double Binary<TLattice, TTraits>::computeEquilibrium(const double& orderparam, const double* velocity, int idx, int k) {

    return Stencil::Weights[idx] * (ChemicalPotential<>::get<TLattice>(k) * mGamma / Stencil::Cs2 + orderparam * CollisionBase<TLattice, Stencil>::computeVelocityFactor(velocity, idx));

}

//FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
//Model is given a "TTraits" class that contains stencil, data, force and boundary information

template<class TLattice>
using DefaultTraitFlowFieldBinary = typename DefaultTrait<TLattice,2> :: template SetBoundary<BounceBack>
                                                                     :: template SetForce<ChemicalForceBinary<Guo,Gradient>>;

template<class TLattice, class TTraits = DefaultTraitFlowFieldBinary<TLattice>>
class FlowFieldBinary : public FlowField<TLattice, TTraits>{ //Inherit from base class to avoid repetition of common
                                                         //calculations

    using Stencil = typename TTraits::Stencil;  
    static constexpr int mNDIM = TLattice::NDIM; 
    
    public:

        inline void setTau1(double val){mTau1 = val;}
        inline void setTau2(double val){mTau2 = val;}

        inline virtual void collide() override; //Collision step

        inline virtual void initialise() override; //Initialisation step

    private:

        inline double computeEquilibrium(const double& density, const double* velocity, const double& order_parameter, const double& chemical_potential, int idx, int k); //Calculate equilibrium in direction idx with a given//density and velocity
                                                                     //at index idx

        double mTau1 = 1;
        double mTau2 = 1;

        enum{x = 0, y = 1, z = 2};

};

template<class TLattice, class TTraits>
inline double FlowFieldBinary<TLattice, TTraits>::computeEquilibrium(const double& density, const double* velocity, const double& order_parameter, const double& chemical_potential, int idx, int k) {

    return density * CollisionBase<TLattice, Stencil>::computeGamma(velocity, idx) + Stencil::Weights[idx] * order_parameter * chemical_potential / Stencil::Cs2; //Equilibrium is density times gamma in this case

}

template<class TLattice, class TTraits>
inline void FlowFieldBinary<TLattice, TTraits>::collide() { //Collision step
       
    #pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k <TLattice::N - TLattice::HaloSize; k++) { //loop over k

        if(!Geometry<TLattice>::isSolid(k)){

            double* old_distribution = this -> mDistribution.getDistributionOldPointer(k);
            double equilibriumsum = 0;
            
            double equilibriums[Stencil::Q];

            for (int idx = 1; idx < Stencil::Q; idx++) {

                equilibriums[idx] = computeEquilibrium(this -> density[k], &(this -> velocity[k * Stencil::D]), OrderParameter<>::get<TLattice>(k), ChemicalPotential<>::get<TLattice>(k), idx, k);
                equilibriumsum += equilibriums[idx];

            }

            equilibriums[0]=this -> density[k]-equilibriumsum;

            this -> collisionQ(equilibriums,old_distribution,InverseTau<>::get<TLattice>(k),k);

        }
        
    }

    FlowField<TLattice, TTraits>::mData.communicateDistribution();

}

template<class TLattice, class TTraits>
inline void FlowFieldBinary<TLattice, TTraits>::initialise() { //Initialise model

    this -> mData.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this -> mt_Forces,mTau1,mTau2);

    #pragma omp parallel for schedule(guided)
    for (int k = TLattice::HaloSize; k <TLattice::N - TLattice::HaloSize; k++) { //loop over k

        double* distribution = this -> mDistribution.getDistributionPointer(k);
        double* old_distribution = this -> mDistribution.getDistributionOldPointer(k);

        Density<>::initialise<TLattice>(1.0, k); //Set density to 1 initially (This will change)
        Velocity<>::initialise<TLattice,TLattice::NDIM>(0.0, k, x);
        if constexpr (TLattice::NDIM >= 2) Velocity<>::initialise<TLattice,TLattice::NDIM>(0.0, k, y);
        if constexpr (TLattice::NDIM == 3) Velocity<>::initialise<TLattice, TLattice::NDIM>(0.0, k, z);

        double equilibriumsum = 0;

        for (int idx = Stencil::Q-1; idx>= 0; idx--) {

            double equilibrium;

            if (idx>0) equilibrium = computeEquilibrium(this -> density[k], &(this -> velocity[k * Stencil::D]), OrderParameter<>::get<TLattice>(k), ChemicalPotential<>::get<TLattice>(k), idx, k);

            else equilibrium = FlowField<TLattice, TTraits>::density[k] - equilibriumsum;

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;
            
            equilibriumsum += equilibrium;
        }
        
    }
    
}
