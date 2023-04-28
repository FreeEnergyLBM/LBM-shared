#ifndef FLOWFIELDBINARY_HEADER
#define FLOWFIELDBINARY_HEADER
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

//Trait class for FlowField Distribution (Navier-Stokes and continuity solver)

template<int NDIM=3>
struct traitFlowFieldBinaryDefault{
    
    using Stencil=std::conditional_t<NDIM==2,D2Q9,D3Q19>; //Here, D refers to the number of cartesian dimensions
                        //and Q refers to the number of discrete velocity directions.
                        //This naming convention is standard in LBM.

    using Boundaries=LatticeTuple<BounceBack>; //This will tell the model which boundaries to apply
    using Forces=LatticeTuple<BodyForce,ChemicalForce>; //This will tell the model which forces to apply
};

template<int ndim=3,class traits=traitFlowFieldBinaryDefault<ndim>>
class FlowFieldBinary:public FlowField<ndim,traits>{ //Inherit from base class to avoid repetition of common
                                                         //calculations
    
    public:
        template<int lx, int ly,int lz>
        FlowFieldBinary(LatticeProperties<lx,ly,lz>& properties):FlowField<ndim,traits>(properties),m_OrderParameter(properties),m_ChemicalPotential(properties){}

        virtual void collide() override; //Collision step

        virtual void initialise() override; //Initialisation step

    private:

        double computeEquilibrium(const double& density,const double* velocity,const double& order_parameter,const double& chemical_potential,const int idx,const int k) const; //Calculate equilibrium in direction idx with a given//density and velocity

        double computeCollisionQ(double& sum,const int k,const double& old,const double& density,
                                            const double* velocity,const double& order_parameter,const double& chemical_potential,const int idx) const; //Calculate collision                                                                             //at index idx


        OrderParameter m_OrderParameter;
        ChemicalPotential m_ChemicalPotential;
        enum{x=0,y=1,z=2};

};

template<int ndim,class traits>
double FlowFieldBinary<ndim,traits>::computeEquilibrium(const double& density,const double* velocity,const double& order_parameter,const double& chemical_potential,const int idx,const int k) const{

    return density*CollisionBase<typename traits::Stencil>::computeGamma(velocity,idx)+traits::Stencil::Weights[idx]*order_parameter*chemical_potential/traits::Stencil::Cs2; //Equilibrium is density
                                                                                        //times gamma in this
                                                                                        //case

}

template<int ndim,class traits>
void FlowFieldBinary<ndim,traits>::collide(){ //Collision step

    
    int QQ=traits::Stencil::Q;
    int QN=traits::Stencil::Q*FlowField<ndim,traits>::N;
    int DN=traits::Stencil::D*FlowField<ndim,traits>::N;
    double* old_distribution=FlowField<ndim,traits>::m_Distribution.getDistributionOldPointer(0);
    double* distribution=FlowField<ndim,traits>::m_Distribution.getDistributionPointer(0);

    #ifdef OMPPARALLEL
    double CollideStartTime=omp_get_wtime();
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=FlowField<ndim,traits>::LY*FlowField<ndim,traits>::LZ*MAXNEIGHBORS;k<FlowField<ndim,traits>::N-MAXNEIGHBORS*FlowField<ndim,traits>::LY*FlowField<ndim,traits>::LZ;k++){ //loop over k

        double sum=0;
        for (int idx=traits::Stencil::Q-1;idx>=0;--idx){ //loop over discrete velocity directions
            //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"
            FlowField<ndim,traits>::m_Distribution.getDistributionPointer(FlowField<ndim,traits>::m_Distribution.streamIndex(k,idx))[idx]=computeCollisionQ(sum,k,old_distribution[k*traits::Stencil::Q+idx],FlowField<ndim,traits>::density[k],&FlowField<ndim,traits>::velocity[k*traits::Stencil::D],m_OrderParameter.getParameter(k),m_ChemicalPotential.getParameter(k),idx);
        }        
        
    }
    #ifdef MPIPARALLEL
    FlowField<ndim,traits>::m_Data.communicateDistribution();
    #endif
    
}

template<int ndim,class traits>
void FlowFieldBinary<ndim,traits>::initialise(){ //Initialise model

    FlowField<ndim,traits>::m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    
    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=FlowField<ndim,traits>::LY*FlowField<ndim,traits>::LZ*MAXNEIGHBORS;k<FlowField<ndim,traits>::N-MAXNEIGHBORS*FlowField<ndim,traits>::LY*FlowField<ndim,traits>::LZ;k++){ //loop over k

        double* distribution=FlowField<ndim,traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution=FlowField<ndim,traits>::m_Distribution.getDistributionOldPointer(k);

        FlowField<ndim,traits>::density[k]=1.0; //Set density to 1 initially (This will change)
        FlowField<ndim,traits>::velocity[k*traits::Stencil::D+x]=0.0; //0 initial velocity
        FlowField<ndim,traits>::velocity[k*traits::Stencil::D+y]=0;
        FlowField<ndim,traits>::velocity[k*traits::Stencil::D+z]=0;
        int sum=0;
        for (int idx=traits::Stencil::Q-1;idx>=0;idx--){

            double equilibrium;
            if (idx>0) equilibrium=computeEquilibrium(FlowField<ndim,traits>::density[k],&FlowField<ndim,traits>::velocity[k*traits::Stencil::D],m_OrderParameter.getParameter(k),m_ChemicalPotential.getParameter(k),idx,k);
            else equilibrium=FlowField<ndim,traits>::density[k]-sum;

            distribution[idx]=equilibrium; //Set distributions to equillibrium
            old_distribution[idx]=equilibrium;        

        }
        
    }
    
}

template<int ndim,class traits>
double FlowFieldBinary<ndim,traits>::computeCollisionQ(double& sum,const int k,const double& old,const double& density,
                                            const double* velocity,const double& order_parameter,const double& chemical_potential,const int idx) const{
                                            //Calculate collision step at a given velocity index at point k
    
    double forcexyz[traits::Stencil::D]; //Temporary array storing force in each cartesian direction

    //Force is the sum of model forces and given forces
    for(int xyz=0;xyz<traits::Stencil::D;xyz++) forcexyz[xyz]=FlowField<ndim,traits>::computeModelForce(xyz,k)+FlowField<ndim,traits>::computeForces(xyz,k);
    
    //Sum of collision + force contributions 
    if (idx>0) {
        double eq=CollisionBase<typename traits::Stencil>::collideSRT(old,computeEquilibrium(density,velocity,order_parameter,chemical_potential,idx,k),FlowField<ndim,traits>::m_InverseTau)
              +CollisionBase<typename traits::Stencil>::forceGuoSRT(forcexyz,velocity,FlowField<ndim,traits>::m_InverseTau,idx);
        sum+=eq;
        
        return eq;
    }
    else return density-sum+CollisionBase<typename traits::Stencil>::forceGuoSRT(forcexyz,velocity,FlowField<ndim,traits>::m_InverseTau,idx);
    
}

#endif