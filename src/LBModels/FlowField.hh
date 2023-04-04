#ifndef FLOWFIELD_HEADER
#define FLOWFIELD_HEADER
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include "../BoundaryModels/Boundaries.hh"
#include "../Forces/Forces.hh"
#include <utility>

//FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
//Model is given a "traits" class that contains stencil, data, force and boundary information

template<class traits>
class FlowField:public CollisionBase<typename traits::Stencil>{ //Inherit from base class to avoid repetition of common
                                                         //calculations
    public:
        //Constructors to construct tuples of forces and boundaries
        FlowField():m_Data(),mt_Forces(*new typename traits::Forces),mt_Boundaries(*new typename traits::Boundaries),m_Distribution(m_Data.getDistributionObject()){
            
        }

        void precompute(); //Perform any necessary computations before collision

        //void collide(); //Collision step

        void boundaries(); //Boundary calculation

        virtual void initialise(); //Initialisation step

        void computeMomenta(); //Momenta (density, velocity) calculation

        const double& getDensity(int k) const; //Return density at lattice point k

        const std::vector<double>& getVelocity() const; //Return vector of velocity

        const std::vector<double>& getDistribution() const; //Return vector of distribution

        template<class traits2>
        friend class FlowFieldBinary;

    private:

       double computeEquilibrium(const double& density,const double* velocity,
                                  const int idx,const int k) const; //Calculate equilibrium in direction idx with a given
                                                        //density and velocity
        //#pragma omp begin declare target
        double computeModelForce(int xyz,int k) const; //Calculate forces specific to the model in direction xyz

        double computeForces(int xyz,int k) const; //Calculate other forces in direction xyz
        //#pragma omp end declare target
        double computeCollisionQ(int k,const double& old,const double& density,
                                 const double* velocity,const int idx) const; //Calculate collision
                                                                                           //at index idx

        double computeDensity(const double* distribution,int k) const; //Calculate density

        double computeVelocity(const double* distribution,const double& density,
                               const int xyz,int k) const; //Calculate velocity
        #pragma omp begin declare target
        static constexpr double m_Tau=1.0; //TEMPORARY relaxation time
        
        static constexpr double m_InverseTau=1.0/m_Tau; //TEMPORARY inverse relaxation time
        #pragma omp end declare target
        Density<double> m_Density; //Density

        Velocity<double,NDIM> m_Velocity; //Velocity

        Distribution_Base<typename traits::Stencil>& m_Distribution;
            //Distributions

        vector<double>& density=m_Density.getParameter(); //Reference to vector of densities

        vector<double>& velocity=m_Velocity.getParameter(); //Reference to vector of velocities

        vector<double>& distribution=m_Distribution.getDistribution(); //Reference to vector of distributions

        enum{x=0,y=1,z=2}; //Indices corresponding to x, y, z directions


        typename traits::Data m_Data; //MOVE THIS TO BASE

        typename traits::Forces& mt_Forces; //MOVE THIS TO BASE
        typename traits::Boundaries& mt_Boundaries;
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
    #pragma omp parallel for schedule( static )
    #endif
    for (int k=LY*LZ*MAXNEIGHBORS;k<N-MAXNEIGHBORS*LY*LZ;k++){ //loop over k

        if constexpr(std::tuple_size<typename traits::Forces>::value!=0){ //Check if there is at least one element
                                                                          //in F
            std::apply([k](auto&... forces){//See Algorithm.hh for explanation of std::apply
                (forces.precompute(k),...);
            }, mt_Forces);
        }
        else;
        //while(m_Geometry.isSolid(k+1)&&k<N-MAXNEIGHBORS*LY*LZ){
        //    k++;
        //}
        //k = m_Data.iterateFluid(k,false); //increment k

    }

    m_Distribution.getDistribution().swap(m_Distribution.getDistributionOld()); //swap old and new distributions
                                                                                //before collision
}
#pragma omp begin declare target
template<class traits>
double FlowField<traits>::computeForces(int xyz,int k) const{ //Return the sum of forces

    if constexpr(std::tuple_size<typename traits::Forces>::value!=0){
        return std::apply([xyz,k](auto&... forces){
                return (forces.compute(xyz,k)+...);
            }, mt_Forces);
    }
    else return 0;

}
#pragma omp end declare target
/*
template<class traits>
void FlowField<traits>::collide(){ //Collision step

    //int k=LY*LZ*MAXNEIGHBORS;
    //k = m_Data.iterateFluid0(k,false);
    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( static )
    #endif
    for (int k=LY*LZ*MAXNEIGHBORS;k<N-MAXNEIGHBORS*LY*LZ;k++){ //loop over k

        double* distribution=m_Distribution.getDistributionPointer(k);
        double* old_distribution=m_Distribution.getDistributionOldPointer(k);

        for (int idx=0;idx<traits::Stencil::Q;idx++){ //loop over discrete velocity directions
            //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"
            m_Distribution.getDistributionPointer(m_Distribution.streamIndex(k,idx))[idx]=computeCollisionQ(k,old_distribution[idx],density[k],&velocity[k*traits::Stencil::D],idx);
        }
        //while(m_Geometry.isSolid(k+1)&&k<N-MAXNEIGHBORS*LY*LZ){
        //    k++;
        //}
        //k = m_Data.iterateFluid(k,false); //increment k
        
    }
    
    m_Data.communicateDistribution();
    
}
*/
template<class traits>
void FlowField<traits>::boundaries(){ //Apply the boundary step

    //int k=0;
    //k = m_Data.iterateSolid0(k,true);
    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( static )
    #endif
    for (int k=0;k<N;k++){ //loop over k
        
        if constexpr(std::tuple_size<typename traits::Boundaries>::value!=0){ //Check if there are any boundary
                                                                              //models
            std::apply([this,k](auto&... boundaries){
                for (int idx=0;idx<traits::Stencil::Q;idx++){
                    if(!m_Geometry.isSolid(m_Distribution.streamIndex(k,idx))){
                        (boundaries.compute(this->m_Distribution,k,idx),...);
                    }
                }
            }, mt_Boundaries);
            
        }
        else;
        //while(!m_Geometry.isSolid(k+1)&&k<N){
        //    k++;
        //}
        //k = m_Data.iterateSolid(k,true); //increment k

    }
    
}

template<class traits>
void FlowField<traits>::initialise(){ //Initialise model

    m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    
    //int k=LY*LZ*MAXNEIGHBORS;
    //k = m_Data.iterateFluid0(k,false);
    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( static )
    #endif
    for (int k=LY*LZ*MAXNEIGHBORS;k<N-MAXNEIGHBORS*LY*LZ;k++){ //loop over k

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
        //while(m_Geometry.isSolid(k+1)&&k<N-MAXNEIGHBORS*LY*LZ){
        //    k++;
        //}
        //k = m_Data.iterateFluid(k,false); //increment k
        
    }
    
}


template<class traits>
void FlowField<traits>::computeMomenta(){ //Calculate Density and Velocity

    //int k=LY*LZ*MAXNEIGHBORS;
    //k = m_Data.iterateFluid0(k,false);
    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( static )
    #endif
    for (int k=LY*LZ*MAXNEIGHBORS;k<N-MAXNEIGHBORS*LY*LZ;k++){ //Loop over k

        double* distribution=m_Distribution.getDistributionPointer(k);

        density[k]=computeDensity(distribution,k); //Calculate density
        velocity[k*traits::Stencil::D+x]=computeVelocity(distribution,density[k],x,k); //Calculate velocities
        velocity[k*traits::Stencil::D+y]=computeVelocity(distribution,density[k],y,k);
        velocity[k*traits::Stencil::D+z]=computeVelocity(distribution,density[k],z,k);
        //while(m_Geometry.isSolid(k+1)&&k<N-MAXNEIGHBORS*LY*LZ){
        //    k++;
        //}
        //k = m_Data.iterateFluid(k,false); //increment k

    }
    
    m_Data.communicate(m_Density);
    //m_Data.communicate(m_Velocity);

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
              +CollisionBase<typename traits::Stencil>::forceSRT(forcexyz,velocity,m_InverseTau,idx);

}


template<class traits>
double FlowField<traits>::computeEquilibrium(const double& density,const double* velocity,const int idx,const int k) const{

    return density*CollisionBase<typename traits::Stencil>::computeGamma(velocity,idx); //Equilibrium is density
                                                                                        //times gamma in this
                                                                                        //case

}
#pragma omp begin declare target
template<class traits>
double FlowField<traits>::computeModelForce(int k,int xyz) const{

    return 0.0; //No model force in this case

}
#pragma omp end declare target
template<class traits>
double FlowField<traits>::computeDensity(const double* distribution,int k) const{ //Density calculation
    //Density is the sum of distributions plus any source/correction terms
    if constexpr(std::tuple_size<typename traits::Forces>::value!=0){
        return CollisionBase<typename traits::Stencil>::computeFirstMoment(distribution)+std::apply([k](auto&... forces){
                return (forces.computeDensitySource(k)+...);
            }, mt_Forces);
    }
    //CHANGE THIS SO FIRST/SECOND MOMENT COMPUTATION IS DONE IN DISTRIBUTION
    else return CollisionBase<typename traits::Stencil>::computeFirstMoment(distribution);

}

template<class traits>
double FlowField<traits>::computeVelocity(const double* distribution,const double& density,
                                          const int xyz,int k) const{ //Velocity calculation in direction xyz
    //Velocity in direction xyz is sum of distribution times the xyz component of the discrete velocity vector
    //in each direction plus any source/correction terms
    if constexpr(std::tuple_size<typename traits::Forces>::value!=0){
        return CollisionBase<typename traits::Stencil>::computeSecondMoment(distribution,xyz)+std::apply([xyz,k](auto&&... forces){
                return (forces.computeVelocitySource(xyz,k)+...);
            }, mt_Forces);
    }
    else return CollisionBase<typename traits::Stencil>::computeSecondMoment(distribution,xyz);

}

#endif