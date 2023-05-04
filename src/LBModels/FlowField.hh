#pragma once
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include "../BoundaryModels/Boundaries.hh"
#include "../Forces/Forces.hh"
#include "../Parallel.hh"
#include <utility>

//FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
//Model is given a "traits" class that contains stencil, data, force and boundary information


//Trait class for FlowField Distribution (Navier-Stokes and continuity solver)
/*
template<int NDIM>
struct traitFlowFieldDefault{

    using Stencil=std::conditional_t<NDIM==2,D2Q9,D3Q19>; //Here, D refers to the number of cartesian dimensions
                        //and Q refers to the number of discrete velocity directions.
                        //This naming convention is standard in LBM.
    using Boundaries=LatticeTuple<BounceBack>; //This will tell the model which boundaries to apply
    using Forces=LatticeTuple<>; //This will tell the model which forces to apply
};
*/

template<class traits>
class FlowField : public CollisionBase<typename traits::Stencil> { //Inherit from base class to avoid repetition of common
                                                         //calculations
    public:

        //Constructors to construct tuples of forces and boundaries
        
        constexpr FlowField(typename traits::Properties& properties)
          : CollisionBase<typename traits::Stencil>(properties),
            m_Density(properties),
            m_Velocity(properties),
            m_Distribution(m_Data.getDistributionObject()),
            N(properties.m_N),
            LX(properties.m_LX),
            LY(properties.m_LY),
            LZ(properties.m_LZ),
            NDIM(properties.m_NDIM),
            HaloSize(properties.m_HaloSize),
            m_Data(properties),
            mt_Forces(properties),
            mt_Boundaries(properties),
            m_Geometry(properties)
        {}

        void precompute(); //Perform any necessary computations before collision

        virtual void collide(); //Collision step

        void boundaries(); //Boundary calculation

        virtual void initialise(); //Initialisation step

        void computeMomenta(); //Momenta (density, velocity) calculation

        const double& getDensity(int k) const; //Return density at lattice point k

        const std::vector<double>& getVelocity() const; //Return vector of velocity

        const std::vector<double>& getDistribution() const; //Return vector of distribution

        template<class>
        friend class FlowFieldBinary;

        template<class force,int inst=0>
        force& getForce() {
            auto forces=get_type<force>(mt_Forces.getTuple());
            return std::get<inst>(forces);
        }

        template<class boundary,int inst=0>
        boundary& getBoundary() {
            auto boundaries=get_type<boundary>(mt_Boundaries.getTuple());
            return std::get<inst>(boundaries);
        }

    private:

       double computeEquilibrium(const double& density,const double* velocity,
                                  const int idx,const int k) const; //Calculate equilibrium in direction idx with a given
                                                        //density and velocity

        double computeModelForce(int xyz,int k) const; //Calculate forces specific to the model in direction xyz

        double computeForces(int xyz,int k) const; //Calculate other forces in direction xyz

        double computeCollisionQ(int k,const double& old,const double& density,
                                 const double* velocity,const int idx) const; //Calculate collision
                                                                                           //at index idx

        double computeDensity(const double* distribution,int k) const; //Calculate density

        double computeVelocity(const double* distribution,const double& density,
                               const int xyz,int k) const; //Calculate velocity

        static constexpr double m_Tau=1.0; //TEMPORARY relaxation time
        
        static constexpr double m_InverseTau=1.0/m_Tau; //TEMPORARY inverse relaxation time

        Density m_Density; //Density

        Velocity m_Velocity; //Velocity

        typename traits::Properties::template DataType<typename traits::Stencil>::DistributionData& m_Distribution;
            //Distributions

        vector<double>& density=m_Density.getParameter(); //Reference to vector of densities

        vector<double>& velocity=m_Velocity.getParameter(); //Reference to vector of velocities

        vector<double>& distribution=m_Distribution.getDistribution(); //Reference to vector of distributions

        enum{x=0,y=1,z=2}; //Indices corresponding to x, y, z directions

        const int& N;
        const int& LX;
        const int& LY;
        const int& LZ;
        const int& NDIM;
        const int& HaloSize;

        typename traits::Properties::template DataType<typename traits::Stencil> m_Data; //MOVE THIS TO BASE

        typename traits::Forces mt_Forces; //MOVE THIS TO BASE
        typename traits::Boundaries mt_Boundaries;
        Geometry m_Geometry;
};


template<class traits>
const double& FlowField<traits>::getDensity(int k) const{

    return density[k]; //Return reference to density at point k

}

template<class traits>
const std::vector<double>& FlowField<traits>::getVelocity() const{

    return velocity; //Return reference to velocity vector

}

template<class traits>
const std::vector<double>& FlowField<traits>::getDistribution() const{

    return distribution; //Return reference to distribution vector

}

template<class traits>
void FlowField<traits>::precompute(){ //Perform necessary calculations before collision

    
    //k = m_Data.iterateFluid0(k,false);
    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=HaloSize;k<N-HaloSize;k++){ //loop over k

        if constexpr(std::tuple_size<typename traits::Forces::getTupleType>::value!=0){ //Check if there is at least one element
                                                                          //in F
            std::apply([k](auto&... forces){//See Algorithm.hh for explanation of std::apply
                (forces.precompute(k),...);
            }, mt_Forces.getTuple());
        }
        else;


    }

    m_Distribution.getDistribution().swap(m_Distribution.getDistributionOld()); //swap old and new distributions
                                                                                //before collision
}

template<class traits>
double FlowField<traits>::computeForces(int xyz,int k) const{ //Return the sum of forces

    if constexpr (std::tuple_size<typename traits::Forces::getTupleType>::value!=0){
        return std::apply([xyz,k](auto&... forces){
                return (forces.compute(xyz,k)+...);
            }, mt_Forces.getTuple());
    }
    else return 0;

}

template<class traits>
void FlowField<traits>::collide(){ //Collision step

    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=HaloSize;k<N-HaloSize;k++){ //loop over k

        double* old_distribution = m_Distribution.getDistributionOldPointer(k);

        for (int idx=0;idx<traits::Stencil::Q;idx++){ //loop over discrete velocity directions
            //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"
            double collision = computeCollisionQ(k, old_distribution[idx], density[k], &velocity[k*traits::Stencil::D], idx);
            m_Distribution.getDistributionPointer(m_Distribution.streamIndex(k,idx))[idx] = collision;
        }
        
    }
    #ifdef MPIPARALLEL
    m_Data.communicateDistribution();
    #endif
}

template<class traits>
void FlowField<traits>::boundaries(){ //Apply the boundary step

    //int k=0;
    //k = m_Data.iterateSolid0(k,true);
    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=0;k<N;k++){ //loop over k
        
        if constexpr(std::tuple_size<typename traits::Boundaries::getTupleType>::value!=0){ //Check if there are any boundary
                                                                              //models
            std::apply([this,k](auto&... boundaries){
                for (int idx=0;idx<traits::Stencil::Q;idx++){
                    if(m_Geometry.isSolid(k)&&!m_Geometry.isSolid(m_Distribution.streamIndex(k,idx))){
                        (boundaries.compute(this->m_Distribution,k,idx),...);
                    }
                }
            }, mt_Boundaries.getTuple());
            
        }
        //while(!m_Geometry.isSolid(k+1)&&k<N){
        //    k++;
        //}
        //k = m_Data.iterateSolid(k,true); //increment k

    }
    
}

template<class traits>
void FlowField<traits>::initialise(){ //Initialise model

    m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    
    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=HaloSize;k<N-HaloSize;k++){ //loop over k

        double* distribution=m_Distribution.getDistributionPointer(k);
        double* old_distribution=m_Distribution.getDistributionOldPointer(k);

        density[k]=1.0; //Set density to 1 initially (This will change)
        velocity[k*traits::Stencil::D+x]=0.0; //0 initial velocity
        velocity[k*traits::Stencil::D+y]=0;
        velocity[k*traits::Stencil::D+z]=0;

        for (int idx=0;idx<traits::Stencil::Q;idx++){

            double equilibrium=computeEquilibrium(density[k],&velocity[k*traits::Stencil::D],idx,k);

            distribution[idx]=equilibrium; //Set distributions to equillibrium
            old_distribution[idx]=equilibrium;        

        }
        
    }
    
}


template<class traits>
void FlowField<traits>::computeMomenta(){ //Calculate Density and Velocity

    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=HaloSize;k<N-HaloSize;k++){ //Loop over k

        double* distribution=m_Distribution.getDistributionPointer(k);

        density[k]=computeDensity(distribution,k); //Calculate density
        velocity[k*traits::Stencil::D+x]=computeVelocity(distribution,density[k],x,k); //Calculate velocities
        velocity[k*traits::Stencil::D+y]=computeVelocity(distribution,density[k],y,k);
        if constexpr (traits::Stencil::D==3)velocity[k*traits::Stencil::D+z]=computeVelocity(distribution,density[k],z,k);
    }
    #ifdef MPIPARALLEL
    m_Data.communicate(m_Density);
    #endif

}

template<class traits>
double FlowField<traits>::computeCollisionQ(const int k,const double& old,const double& density,
                                            const double* velocity,const int idx) const{
                                            //Calculate collision step at a given velocity index at point k

    double forcexyz[traits::Stencil::D]; //Temporary array storing force in each cartesian direction

    //Force is the sum of model forces and given forces
    for(int xyz=0;xyz<traits::Stencil::D;xyz++) forcexyz[xyz]=computeModelForce(xyz,k)+computeForces(xyz,k);
    
    //Sum of collision + force contributions
    return CollisionBase<typename traits::Stencil>::collideSRT(old,computeEquilibrium(density,velocity,idx,k),m_InverseTau)
              +CollisionBase<typename traits::Stencil>::forceGuoSRT(forcexyz,velocity,m_InverseTau,idx);

}


template<class traits>
double FlowField<traits>::computeEquilibrium(const double& density,const double* velocity,const int idx,const int k) const{

    return density*CollisionBase<typename traits::Stencil>::computeGamma(velocity,idx); //Equilibrium is density
                                                                                        //times gamma in this
                                                                                        //case

}

template<class traits>
double FlowField<traits>::computeModelForce(int k,int xyz) const{

    return 0.0; //No model force in this case

}

template<class traits>
double FlowField<traits>::computeDensity(const double* distribution,int k) const{ //Density calculation
    //Density is the sum of distributions plus any source/correction terms
    if constexpr(std::tuple_size<typename traits::Forces::getTupleType>::value!=0){
        return CollisionBase<typename traits::Stencil>::computeZerothMoment(distribution)+std::apply([k](auto&... forces){
                return (forces.computeDensitySource(k)+...);
            }, mt_Forces.getTuple());
    }
    //CHANGE THIS SO FIRST/SECOND MOMENT COMPUTATION IS DONE IN DISTRIBUTION
    else return CollisionBase<typename traits::Stencil>::computeZerothMoment(distribution);

}

template<class traits>
double FlowField<traits>::computeVelocity(const double* distribution,const double& density,
                                          const int xyz,int k) const{ //Velocity calculation in direction xyz
    //Velocity in direction xyz is sum of distribution times the xyz component of the discrete velocity vector
    //in each direction plus any source/correction terms
    if constexpr(std::tuple_size<typename traits::Forces::getTupleType>::value!=0){
        return CollisionBase<typename traits::Stencil>::computeFirstMoment(distribution,xyz)+std::apply([xyz,k](auto&&... forces){
                return (forces.computeVelocitySource(xyz,k)+...);
            }, mt_Forces.getTuple());
    }
    else return CollisionBase<typename traits::Stencil>::computeFirstMoment(distribution,xyz);

}

template<class properties>
struct DefaultTrait<properties,FlowField>{
    using Stencil=std::conditional_t<properties::m_NDIM==2,D2Q9,D3Q19>; //Here, D refers to the number of cartesian dimensions
                        //and Q refers to the number of discrete velocity directions.
                        //This naming convention is standard in LBM.
    using Boundaries=LatticeTuple<BounceBack>; //This will tell the model which boundaries to apply
    using Forces=LatticeTuple<>; //This will tell the model which forces to apply
    using Properties=properties;
};
