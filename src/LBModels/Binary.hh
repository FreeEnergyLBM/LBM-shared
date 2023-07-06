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
//Model is given a "traits" class that contains stencil, data, force and boundary information

template<class lattice>
using DefaultTraitBinary = typename DefaultTrait<lattice> :: template SetBoundary<BounceBack> 
                                                          :: template SetPreProcessor<ChemicalPotentialCalculatorBinary,CubicWetting> 
                                                          :: template SetPostProcessor<GradientsMultiStencil<OrderParameter<>,CentralXYZ,LaplacianCentral>>;


template<class lattice, class traits = DefaultTraitBinary<lattice>>
class Binary: public CollisionBase<lattice, typename traits::Stencil>, public ModelBase<lattice, traits> { //Inherit from base class to avoid repetition of common
                                                      //calculations
                                                      
    public:

        inline void setTau1(double val){m_Tau1=val;}
        inline void setTau2(double val){m_Tau2=val;}

        inline void collide() override; //Collision step

        inline void initialise() override; //Initialisation step

        inline void computeMomenta() override; //Momenta (density, velocity) calculation

        inline const double& getDensity(const int k) const; //Return density at lattice point k

        inline const std::vector<double>& getVelocity() const; //Return vector of velocity

    private:

        inline double computeEquilibrium(const double& orderparam, const double* velocity,
                                   const int idx, const int k) const; //Calculate equilibrium in direction idx with a given
                                                        //density and velocity

        inline double computeModelForce(int xyz, int k) const; //Calculate forces specific to the model in direction xyz

        inline double computeCollisionQ(const double* forcemethods,const double* equilibriums, int k, const double* old, const int idx) const; //Calculate collision
                                                                                           //at index idx

        inline double computeOrderParameter(const double* distribution, const int k) const; //Calculate the order parameter
                                                                              //corresponding to the relative
                                                                              //concentrations of each phase

        static constexpr double m_Tau = 1.0; //TEMPORARY relaxation time
        static constexpr double m_InverseTau = 1.0 / m_Tau; //TEMPORARY inverse relaxation time

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions
        
        double m_Gamma = 1;

        double m_Tau1=1;
        double m_Tau2=1;

        std::vector<double>& orderparameter = OrderParameter<>::get<typename traits::Lattice>(); //Reference to vector of order parameters
        std::vector<double>& velocity = Velocity<>::get<typename traits::Lattice,traits::Lattice::m_NDIM>(); //Reference to vector of velocities
        std::vector<double>& itau = InverseTau<>::get<typename traits::Lattice>(); //Reference to vector of velocities

};

template<class lattice, class traits>
inline const double& Binary<lattice, traits>::getDensity(const int k) const { //This needs to be renamed

    return orderparameter[k]; //Return reference to order parameter at point k

}

template<class lattice, class traits>
inline const std::vector<double>& Binary<lattice, traits>::getVelocity() const {

    return velocity; //Return reference to velocity vector

}

template<class lattice, class traits>
inline void Binary<lattice, traits>::collide() {

    #pragma omp for schedule(guided)
    for (int k = lattice::m_HaloSize; k < lattice::m_N - lattice::m_HaloSize; k++){ //loop over k
        
        if(!ModelBase<lattice, traits>::m_Geometry.isSolid(k)){
            
            auto forcemethods=ModelBase<lattice,traits>::getForceCalculator(ModelBase<lattice,traits>::mt_Forces,k);

            double* old_distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionOldPointer(k);

            double equilibriumsum = 0;

            double equilibriums[traits::Stencil::Q] = {};
            double forces[traits::Stencil::Q] = {};

            for (int idx = 1; idx <traits::Stencil::Q; idx++) {
                equilibriums[idx]=computeEquilibrium(orderparameter[k], &velocity[k * traits::Stencil::D], idx, k);
                equilibriumsum+=equilibriums[idx];
                if constexpr(std::tuple_size<typename std::remove_reference<decltype(*forcemethods)>::type>::value != 0){
                    forces[idx]=std::apply([idx, k](auto&... forcetype){
                                return (forcetype.template compute<traits>(idx, k) + ...);
                            }, *forcemethods);
                }
            }

            equilibriums[0]=orderparameter[k]-equilibriumsum;

            if constexpr(std::tuple_size<typename std::remove_reference<decltype(*forcemethods)>::type>::value != 0){
                forces[0]=std::apply([k](auto&... forcetype){
                            return (forcetype.template compute<traits>(0, k) + ...);
                        }, *forcemethods);
            }

            this->collisionQ(forces,equilibriums,old_distribution,m_InverseTau,k);

        }
        
    }

    ModelBase<lattice, traits>::m_Data.communicateDistribution();

}

template<class lattice, class traits>
inline void Binary<lattice, traits>::initialise() { //Initialise model

    ModelBase<lattice, traits>::m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    traits::template CollisionModel<typename traits::Stencil>::template initialise<typename traits::Lattice>(m_Tau1,m_Tau2);
    
    #pragma omp parallel for schedule(guided) 
    for (int k = lattice::m_HaloSize; k<lattice::m_N - lattice::m_HaloSize; k++) { //loop over k

        double* distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionOldPointer(k);

        ChemicalPotential<>::initialise<typename traits::Lattice>(0,k);

        OrderParameter<>::initialise<typename traits::Lattice>(1.0,k);
        
        InverseTau<>::initialise<typename traits::Lattice>(1.0/(0.5*(1.0+OrderParameter<>::get<typename traits::Lattice>(k))*(m_Tau1-m_Tau2)+m_Tau2),k);

        double equilibriumsum = 0;

        for (int idx = traits::Stencil::Q - 1; idx>= 0; idx--) {

            double equilibrium;

            if (idx> 0) equilibrium = computeEquilibrium(orderparameter[k], &velocity[k * traits::Stencil::D], idx, k);
            else equilibrium = orderparameter[k] - equilibriumsum;

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;        

            equilibriumsum+=equilibrium;

        }
        
    }

    #ifdef MPIPARALLEL
    #pragma omp master
    {
    ModelBase<lattice, traits>::m_Data.communicate(SolidLabels<>::getInstance<typename traits::Lattice>());
    }
    #endif
    
}


template<class lattice, class traits>
inline void Binary<lattice, traits>::computeMomenta() { //Calculate order parameter

    #pragma omp for schedule(guided)
    for (int k = lattice::m_HaloSize; k <lattice::m_N - lattice::m_HaloSize; k++) { //Loop over k

        if(!Geometry<lattice>::isSolid(k)){

            double* distribution = ModelBase<lattice, traits>::m_Distribution.getDistributionPointer(k);

            orderparameter[k] = computeOrderParameter(distribution, k);

            itau[k]=1.0/(0.5*(1.0+orderparameter[k])*(m_Tau1)-0.5*(-1.0+orderparameter[k])*m_Tau2);
        
        }

    }

    ModelBase<lattice, traits>::m_Data.communicate(OrderParameter<>::getInstance<typename traits::Lattice>());
}

template<class lattice, class traits>
inline double Binary<lattice, traits>::computeEquilibrium(const double& orderparam, const double* velocity, const int idx, const int k) const {

    return traits::Stencil::Weights[idx] * (ChemicalPotential<>::get<typename traits::Lattice>(k) * m_Gamma / traits::Stencil::Cs2 + orderparam * CollisionBase<lattice, typename traits::Stencil>::computeVelocityFactor(velocity, idx));

}

template<class lattice, class traits>
inline double Binary<lattice, traits>::computeOrderParameter(const double* distribution, const int k) const {//Order parameter calculation
    //Order parameter is the sum of distributions plus any source/correction terms

    if constexpr(std::tuple_size<typename traits::Forces>::value != 0){

        return CollisionBase<lattice, typename traits::Stencil>::computeZerothMoment(distribution)
        + std::apply([k](auto&... forces) {

            return (forces.template computeDensitySource<traits>(k) + ...);

        }, ModelBase<lattice, traits>::mt_Forces);

    }
    else return CollisionBase<lattice, typename traits::Stencil>::computeZerothMoment(distribution);

}

//FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
//Model is given a "traits" class that contains stencil, data, force and boundary information

template<class lattice>
using DefaultTraitFlowFieldBinary = typename DefaultTrait<lattice> :: template SetBoundary<BounceBack>
                                                                   :: template SetForce<ChemicalForce<Guo,Gradient>>;

template<class lattice, class traits = DefaultTraitFlowFieldBinary<lattice>>
class FlowFieldBinary : public FlowField<lattice, traits>{ //Inherit from base class to avoid repetition of common
                                                         //calculations
    
    public:

        inline void setTau1(double val){m_Tau1=val;}
        inline void setTau2(double val){m_Tau2=val;}

        inline virtual void collide() override; //Collision step

        inline virtual void initialise() override; //Initialisation step

    private:

        inline double computeEquilibrium(const double& density, const double* velocity, const double& order_parameter, const double& chemical_potential, const int idx, const int k) const; //Calculate equilibrium in direction idx with a given//density and velocity
                                                                     //at index idx

        double m_Tau1=1;
        double m_Tau2=1;

        enum{x=0,y=1,z=2};

};

template<class lattice, class traits>
inline double FlowFieldBinary<lattice, traits>::computeEquilibrium(const double& density, const double* velocity, const double& order_parameter, const double& chemical_potential, const int idx, const int k) const {

    return density * CollisionBase<lattice, typename traits::Stencil>::computeGamma(velocity, idx) + traits::Stencil::Weights[idx] * order_parameter * chemical_potential / traits::Stencil::Cs2; //Equilibrium is density times gamma in this case

}

template<class lattice, class traits>
inline void FlowFieldBinary<lattice, traits>::collide() { //Collision step
       
    #pragma omp for schedule(guided)
    for (int k = lattice::m_HaloSize; k <lattice::m_N - lattice::m_HaloSize; k++) { //loop over k

        if(!Geometry<typename traits::Lattice>::isSolid(k)){

            auto forcemethods=ModelBase<lattice,traits>::getForceCalculator(ModelBase<lattice,traits>::mt_Forces,k);

            //int QQ = traits::Stencil::Q;
            double* old_distribution = FlowField<lattice, traits>::m_Distribution.getDistributionOldPointer(k);
            double equilibriumsum = 0;
            
            double equilibriums[traits::Stencil::Q];
            double forces[traits::Stencil::Q];

            for (int idx = 1; idx <traits::Stencil::Q; idx++) {
                equilibriums[idx]=computeEquilibrium(FlowField<lattice, traits>::density[k],&FlowField<lattice, traits>::velocity[k * traits::Stencil::D], OrderParameter<>::get<typename traits::Lattice>(k), ChemicalPotential<>::get<typename traits::Lattice>(k), idx, k);
                equilibriumsum+=equilibriums[idx];
                if constexpr(std::tuple_size<typename std::remove_reference<decltype(*forcemethods)>::type>::value != 0){
                    forces[idx]=std::apply([idx, k](auto&... forcetype){
                                return (forcetype.template compute<traits>(idx, k) + ...);
                            }, *forcemethods);
                }
            }

            if constexpr(std::tuple_size<typename std::remove_reference<decltype(*forcemethods)>::type>::value != 0){
                forces[0]=std::apply([k](auto&... forcetype){
                            return (forcetype.template compute<traits>(0, k) + ...);
                        }, *forcemethods);
            }

            equilibriums[0]=FlowField<lattice, traits>::density[k]-equilibriumsum;

            this->collisionQ(forces,equilibriums,old_distribution,InverseTau<>::get<typename traits::Lattice>(k),k);

        }
        
    }

    FlowField<lattice, traits>::m_Data.communicateDistribution();
}

template<class lattice, class traits>
inline void FlowFieldBinary<lattice, traits>::initialise() { //Initialise model

    FlowField<lattice, traits>::m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    traits::template CollisionModel<typename traits::Stencil>::template initialise<typename traits::Lattice>(m_Tau1,m_Tau2);

    #pragma omp parallel for schedule(guided)
    for (int k = lattice::m_HaloSize; k <lattice::m_N - lattice::m_HaloSize; k++) { //loop over k

        double* distribution = FlowField<lattice, traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution = FlowField<lattice, traits>::m_Distribution.getDistributionOldPointer(k);

        Density<>::initialise<typename traits::Lattice>(1.0,k); //Set density to 1 initially (This will change)
        Velocity<>::initialise<typename traits::Lattice,traits::Lattice::m_NDIM>(0.0,k,x);
        Velocity<>::initialise<typename traits::Lattice,traits::Lattice::m_NDIM>(0.0,k,y);
        if constexpr (lattice::m_NDIM==3) Velocity<>::initialise<typename traits::Lattice,traits::Lattice::m_NDIM>(0.0,k,z);

        double equilibriumsum = 0;

        for (int idx = traits::Stencil::Q-1; idx>= 0; idx--) {

            double equilibrium;

            if (idx>0) equilibrium = computeEquilibrium(FlowField<lattice, traits>::density[k],&FlowField<lattice, traits>::velocity[k * traits::Stencil::D], OrderParameter<>::get<typename traits::Lattice>(k), ChemicalPotential<>::get<typename traits::Lattice>(k), idx, k);

            else equilibrium = FlowField<lattice, traits>::density[k] - equilibriumsum;

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;
            
            equilibriumsum+=equilibrium;
        }
        
    }
    
}
