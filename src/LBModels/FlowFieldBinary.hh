#pragma once
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include "../BoundaryModels/Boundaries.hh"
#include "../Forces/Forces.hh"
#include "../Parallel.hh"
#include "../GradientStencils/GradientStencils.hh"
#include <utility>
#include <array>
#include <omp.h>

//FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
//Model is given a "traits" class that contains stencil, data, force and boundary information

template<class lattice>
struct DefaultTraitFlowFieldBinary{

    using Stencil = std::conditional_t<std::remove_reference<lattice>::type::m_NDIM == 2, D2Q9, D3Q19>; //Here, D refers to the number of cartesian dimensions

    using Boundaries = std::tuple<BounceBack<lattice>>;

    using Forces = std::tuple<BodyForce<lattice>,ChemicalForce<lattice>>;

};


template<class lattice, class traits = DefaultTraitFlowFieldBinary<lattice>>
class FlowFieldBinary : public FlowField<lattice, traits>{ //Inherit from base class to avoid repetition of common
                                                         //calculations
    
    public:

        inline virtual void collide() override; //Collision step

        inline virtual void initialise() override; //Initialisation step

    private:

        inline double computeEquilibrium(const double& density, const double* velocity, const double& order_parameter, const double& chemical_potential, const int idx, const int k) const; //Calculate equilibrium in direction idx with a given//density and velocity

        inline double computeCollisionQ(double& sum, const int k, const double& old, const double& density,
                                  const double* velocity, const double& order_parameter, const double& chemical_potential, const int idx) const; //Calculate collision                                                                             //at index idx


        OrderParameter<lattice> m_OrderParameter;
        ChemicalPotential<lattice> m_ChemicalPotential;
        InverseTau<lattice> m_InvTau; 

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
        int QQ = traits::Stencil::Q;
        double* old_distribution = FlowField<lattice, traits>::m_Distribution.getDistributionOldPointer(0);
        double equilibriumsum = 0;

        for (int idx = QQ-1; idx>= 0; --idx) { //loop over discrete velocity directions
            //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"

            double& density = FlowField<lattice, traits>::density[k];
            double *velocity = &FlowField<lattice, traits>::velocity[k * traits::Stencil::D];
            double& orderParameter = m_OrderParameter.getParameter(k);
            double& chemicalPotential = m_ChemicalPotential.getParameter(k);

            double collision = computeCollisionQ(equilibriumsum, k, old_distribution[k*QQ+idx], density, velocity, orderParameter, chemicalPotential, idx);
            FlowField<lattice, traits>::m_Distribution.getDistributionPointer(FlowField<lattice, traits>::m_Distribution.streamIndex(k, idx))[idx] = collision;

        }        
        
    }

    #ifdef MPIPARALLEL
    #pragma omp master
    {
    FlowField<lattice, traits>::m_Data.communicateDistribution();
    }
    
    #endif
    
}

template<class lattice, class traits>
inline void FlowFieldBinary<lattice, traits>::initialise() { //Initialise model

    FlowField<lattice, traits>::m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)

    #pragma omp parallel for schedule(guided)
    for (int k = 0; k<lattice::m_N; k++) { //loop over k

        double* distribution = FlowField<lattice, traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution = FlowField<lattice, traits>::m_Distribution.getDistributionOldPointer(k);

        FlowField<lattice, traits>::m_Density.initialise(1.0,k); //Set density to 1 initially (This will change)
        FlowField<lattice, traits>::m_Velocity.initialise(0.0,k,x);
        FlowField<lattice, traits>::m_Velocity.initialise(0.0,k,x);
        FlowField<lattice, traits>::m_Velocity.initialise(0.0,k,x);

        int equilibriumsum = 0;
        for (int idx = traits::Stencil::Q-1; idx>= 0; idx--) {

            double equilibrium;

            if (idx>0) equilibrium = computeEquilibrium(FlowField<lattice, traits>::density[k],&FlowField<lattice, traits>::velocity[k * traits::Stencil::D], m_OrderParameter.getParameter(k), m_ChemicalPotential.getParameter(k), idx, k);

            else equilibrium = FlowField<lattice, traits>::density[k] - equilibriumsum;

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;        

        }
        
    }
    
}

template<class lattice, class traits>
inline double FlowFieldBinary<lattice, traits>::computeCollisionQ(double& equilibriumsum, const int k, const double& old, const double& density,
                                                     const double* velocity, const double& order_parameter, const double& chemical_potential, const int idx) const {
                                            //Calculate collision step at a given velocity index at point k
    
    double forcexyz[traits::Stencil::D]; //Temporary array storing force in each cartesian direction

    //Force is the sum of model forces and given forces
    for(int xyz = 0; xyz <traits::Stencil::D; xyz++) forcexyz[xyz] = FlowField<lattice, traits>::computeModelForce(xyz, k)+FlowField<lattice, traits>::computeForces(xyz, k);
    
    //Sum of collision + force contributions 
    if (idx>0) {

        double eq = CollisionBase<lattice, typename traits::Stencil>::collideSRT(old, computeEquilibrium(density ,velocity ,order_parameter ,chemical_potential ,idx ,k), m_InvTau.getParameter(k))
              + CollisionBase<lattice, typename traits::Stencil>::forceGuoSRT(forcexyz, velocity, m_InvTau.getParameter(k), idx);

        equilibriumsum += eq;
        
        return eq;

    }
    else return density-equilibriumsum+CollisionBase<lattice, typename traits::Stencil>::forceGuoSRT(forcexyz, velocity, m_InvTau.getParameter(k), idx);
    
}
