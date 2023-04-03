#ifndef FLOWFIELDBINARY_HEADER
#define FLOWFIELDBINARY_HEADER
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include "../BoundaryModels/Boundaries.hh"
#include "../Forces/Forces.hh"
#include <utility>
#include <array>

//FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
//Model is given a "traits" class that contains stencil, data, force and boundary information

template<class traits>
class FlowFieldBinary:public FlowField<traits>{ //Inherit from base class to avoid repetition of common
                                                         //calculations
    
    public:
        FlowFieldBinary(typename traits::Forces& forces,typename traits::Boundaries& boundaries):FlowField<traits>(forces,boundaries){}
        FlowFieldBinary(typename traits::Forces& forces):FlowField<traits>(forces){}
        FlowFieldBinary(typename traits::Boundaries& boundaries):FlowField<traits>(boundaries){}
        FlowFieldBinary():FlowField<traits>(){}

        virtual void collide() override; //Collision step

        virtual void initialise() override; //Initialisation step

    private:
        virtual double computeEquilibrium(const double& density,const double* velocity,
                                  const int idx,const int k) const override; //Calculate equilibrium in direction idx with a given
                                                        //density and velocity

        double computeCollisionQ(double& sum, int k,const double& old,const double& density,
                                 const double* velocity,const int idx) const; //Calculate collision
                                                                                           //at index idx
        

        OrderParameter<double> m_OrderParameter;
        ChemicalPotential<double> m_ChemicalPotential;
        enum{x=0,y=1,z=2};

};

template<class traits>
double FlowFieldBinary<traits>::computeEquilibrium(const double& density,const double* velocity,const int idx,const int k) const{

    return density*CollisionBase<typename traits::Stencil>::computeGamma(velocity,idx)+traits::Stencil::Weights[idx]*m_OrderParameter.getParameter(k)*m_ChemicalPotential.getParameter(k)/traits::Stencil::Cs2; //Equilibrium is density
                                                                                        //times gamma in this
                                                                                        //case

}

template<class traits>
void FlowFieldBinary<traits>::collide(){ //Collision step

    //int k=LY*LZ*MAXNEIGHBORS;
    //k = FlowField<traits>::m_Data.iterateFluid0(k,false);
    #ifdef OMPPARALLEL
    #pragma omp for
    #endif
    for (int k=LY*LZ*MAXNEIGHBORS;k<N-MAXNEIGHBORS*LY*LZ;k++){ //loop over k

        double* distribution=FlowField<traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution=FlowField<traits>::m_Distribution.getDistributionOldPointer(k);
        double sum=0;
        for (int idx=traits::Stencil::Q-1;idx>=0;--idx){ //loop over discrete velocity directions
            //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"
            FlowField<traits>::m_Distribution.getDistributionPointer(FlowField<traits>::m_Distribution.streamIndex(k,idx))[idx]=computeCollisionQ(sum,k,old_distribution[idx],FlowField<traits>::density[k],&FlowField<traits>::velocity[k*traits::Stencil::D],idx);
        }
        while(FlowField<traits>::m_Geometry.isSolid(k+1)&&k<N-MAXNEIGHBORS*LY*LZ){
            k++;
        }
        //k = FlowField<traits>::m_Data.iterateFluid(k,false); //increment k
        
    }
    
    FlowField<traits>::m_Data.communicateDistribution();
    
}

template<class traits>
void FlowFieldBinary<traits>::initialise(){ //Initialise model

    FlowField<traits>::m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    
    //int k=LY*LZ*MAXNEIGHBORS;
    //k = FlowField<traits>::m_Data.iterateFluid0(k,false);
    #ifdef OMPPARALLEL
    #pragma omp for
    #endif
    for (int k=LY*LZ*MAXNEIGHBORS;k<N-MAXNEIGHBORS*LY*LZ;k++){ //loop over k

        double* distribution=FlowField<traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution=FlowField<traits>::m_Distribution.getDistributionOldPointer(k);

        FlowField<traits>::density[k]=1.0; //Set density to 1 initially (This will change)
        FlowField<traits>::velocity[k*traits::Stencil::D+x]=0.0; //0 initial velocity
        FlowField<traits>::velocity[k*traits::Stencil::D+y]=0;
        FlowField<traits>::velocity[k*traits::Stencil::D+z]=0;
        int sum=0;
        for (int idx=traits::Stencil::Q-1;idx>=0;idx--){

            double equilibrium;
            if (idx>0) equilibrium=computeEquilibrium(FlowField<traits>::density[k],&FlowField<traits>::velocity[k*traits::Stencil::D],idx,k);
            else equilibrium=FlowField<traits>::density[k]-sum;

            distribution[idx]=equilibrium; //Set distributions to equillibrium
            old_distribution[idx]=equilibrium;        

        }
        while(FlowField<traits>::m_Geometry.isSolid(k+1)&&k<N-MAXNEIGHBORS*LY*LZ){
            k++;
        }
        //k = FlowField<traits>::m_Data.iterateFluid(k,false); //increment k
        
    }
    
}

template<class traits>
double FlowFieldBinary<traits>::computeCollisionQ(double& sum,const int k,const double& old,const double& density,
                                            const double* velocity,const int idx) const{
                                            //Calculate collision step at a given velocity index at point k

    std::array<double,traits::Stencil::D> forcexyz; //Temporary array storing force in each cartesian direction

    //Force is the sum of model forces and given forces
    for(int xyz=0;xyz<traits::Stencil::D;xyz++) forcexyz[xyz]=FlowField<traits>::computeModelForce(xyz,k)+FlowField<traits>::computeForces(xyz,k);
    
    //Sum of collision + force contributions
    if (idx>0) {
        double eq=CollisionBase<typename traits::Stencil>::collideSRT(old,computeEquilibrium(density,velocity,idx,k),FlowField<traits>::m_InverseTau)
              +CollisionBase<typename traits::Stencil>::forceSRT(forcexyz,velocity,FlowField<traits>::m_InverseTau,idx);
        sum+=eq;
        
        return eq;
    }
    else return FlowField<traits>::m_Density.getParameter(k)-sum+CollisionBase<typename traits::Stencil>::forceSRT(forcexyz,velocity,FlowField<traits>::m_InverseTau,idx);

}



#endif