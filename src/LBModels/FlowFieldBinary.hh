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

template<class properties>
struct DefaultTraitFlowFieldBinary{

    using Stencil = std::conditional_t<std::remove_reference<properties>::type::m_NDIM == 2, D2Q9, D3Q19>; //Here, D refers to the number of cartesian dimensions

    using Boundaries = std::tuple<BounceBack>;

    using Forces = std::tuple<BodyForce,ChemicalForce>;

    using Properties = properties;

};


template<class traits = DefaultTraitFlowFieldBinary<decltype(GETPROPERTIES())>>
class FlowFieldBinary : public FlowField<traits>{ //Inherit from base class to avoid repetition of common
                                                         //calculations
    
    public:

        inline virtual void collide() override; //Collision step

        inline virtual void initialise() override; //Initialisation step

    private:

        inline double computeEquilibrium(const double& density, const double* velocity, const double& order_parameter, const double& chemical_potential, const int idx, const int k) const; //Calculate equilibrium in direction idx with a given//density and velocity

        inline double computeCollisionQ(double& sum, const int k, const double& old, const double& density,
                                  const double* velocity, const double& order_parameter, const double& chemical_potential, const int idx) const; //Calculate collision                                                                             //at index idx


        OrderParameter m_OrderParameter;
        ChemicalPotential m_ChemicalPotential;

        enum{x=0,y=1,z=2};

};

template<class traits>
inline double FlowFieldBinary<traits>::computeEquilibrium(const double& density, const double* velocity, const double& order_parameter, const double& chemical_potential, const int idx, const int k) const {

    return density * CollisionBase<typename traits::Stencil>::computeGamma(velocity, idx) + traits::Stencil::Weights[idx] * order_parameter * chemical_potential / traits::Stencil::Cs2; //Equilibrium is density times gamma in this case

}

template<class traits>
inline void FlowFieldBinary<traits>::collide() { //Collision step
    
    int QQ = traits::Stencil::Q;
    double* old_distribution = FlowField<traits>::m_Distribution.getDistributionOldPointer(0);

    #pragma omp parallel for schedule(guided)
    for (int k = GETPROPERTIES().m_HaloSize; k <GETPROPERTIES().m_N - GETPROPERTIES().m_HaloSize; k++) { //loop over k

        double equilibriumsum = 0;

        for (int idx = QQ-1; idx>= 0; --idx) { //loop over discrete velocity directions
            //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"

            double& density = FlowField<traits>::density[k];
            double *velocity = &FlowField<traits>::velocity[k * traits::Stencil::D];
            double& orderParameter = m_OrderParameter.getParameter(k);
            double& chemicalPotential = m_ChemicalPotential.getParameter(k);

            double collision = computeCollisionQ(equilibriumsum, k, old_distribution[k*QQ+idx], density, velocity, orderParameter, chemicalPotential, idx);
            FlowField<traits>::m_Distribution.getDistributionPointer(FlowField<traits>::m_Distribution.streamIndex(k, idx))[idx] = collision;

        }        
        
    }

    #ifdef MPIPARALLEL
    FlowField<traits>::m_Data.communicateDistribution();
    #endif
    
}

template<class traits>
inline void FlowFieldBinary<traits>::initialise() { //Initialise model

    FlowField<traits>::m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)

    #pragma omp parallel for schedule(guided)
    for (int k = GETPROPERTIES().m_HaloSize; k <GETPROPERTIES().m_N - GETPROPERTIES().m_HaloSize; k++) { //loop over k

        double* distribution = FlowField<traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution = FlowField<traits>::m_Distribution.getDistributionOldPointer(k);

        FlowField<traits>::density[k] = 1.0; //Set density to 1 initially (This will change)
        FlowField<traits>::velocity[k * traits::Stencil::D + x] = 0.0; //0 initial velocity
        FlowField<traits>::velocity[k * traits::Stencil::D + y] = 0;
        FlowField<traits>::velocity[k * traits::Stencil::D + z] = 0;
        int equilibriumsum = 0;
        for (int idx = traits::Stencil::Q-1; idx>= 0; idx--) {

            double equilibrium;

            if (idx>0) equilibrium = computeEquilibrium(FlowField<traits>::density[k],&FlowField<traits>::velocity[k * traits::Stencil::D], m_OrderParameter.getParameter(k), m_ChemicalPotential.getParameter(k), idx, k);

            else equilibrium = FlowField<traits>::density[k] - equilibriumsum;

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;        

        }
        
    }
    
}

template<class traits>
inline double FlowFieldBinary<traits>::computeCollisionQ(double& equilibriumsum, const int k, const double& old, const double& density,
                                                     const double* velocity, const double& order_parameter, const double& chemical_potential, const int idx) const {
                                            //Calculate collision step at a given velocity index at point k
    
    double forcexyz[traits::Stencil::D]; //Temporary array storing force in each cartesian direction

    //Force is the sum of model forces and given forces
    for(int xyz = 0; xyz <traits::Stencil::D; xyz++) forcexyz[xyz] = FlowField<traits>::computeModelForce(xyz, k)+FlowField<traits>::computeForces(xyz, k);
    
    //Sum of collision + force contributions 
    if (idx>0) {

        double eq = CollisionBase<typename traits::Stencil>::collideSRT(old, computeEquilibrium(density ,velocity ,order_parameter ,chemical_potential ,idx ,k), FlowField<traits>::m_InverseTau)
              + CollisionBase<typename traits::Stencil>::forceGuoSRT(forcexyz, velocity, FlowField<traits>::m_InverseTau, idx);

        equilibriumsum += eq;
        
        return eq;

    }
    else return density-equilibriumsum+CollisionBase<typename traits::Stencil>::forceGuoSRT(forcexyz, velocity, FlowField<traits>::m_InverseTau, idx);
    
}
