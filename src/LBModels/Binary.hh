#ifndef BINARY_HEADER
#define BINARY_HEADER
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include "../Parallel.hh"
#include <utility>

//Binary.hh: Contains the details of the LBM model to solve an equation for phase separation. Each
//Model is given a "traits" class that contains stencil, data, force and boundary information


//Trait class for PhaseField Distribution (Calculates the interface between components)
template<int NDIM=3>
struct traitBinaryDefault{

    using Stencil=std::conditional_t<NDIM==2,D2Q9,D3Q19>; //Here, D refers to the number of cartesian dimensions
    using Boundaries=LatticeTuple<BounceBack>;
    using Forces=LatticeTuple<OrderParameterGradients<CentralXYZ<Stencil,ParallelType>>>;
};

template<int ndim=3,class traits=traitBinaryDefault<ndim>>
class Binary:CollisionBase<typename traits::Stencil>{ //Inherit from base class to avoid repetition of common
                                                      //calculations
    public:
        //Constructors to construct tuples of forces and boundaries
        template<int lx, int ly,int lz>
        Binary(LatticeProperties<lx,ly,lz>& properties):LX(properties.m_LX),LY(properties.m_LY),LZ(properties.m_LZ),N(properties.m_N),m_Data(properties),mt_Forces(properties),mt_Boundaries(properties),m_ChemicalPotential(properties),m_GradOrderParameter(properties),m_LaplacianOrderParameter(properties),m_OrderParameter(properties),m_Velocity(properties),CollisionBase<typename traits::Stencil>(properties),m_Geometry(properties){
            
        }

        void precompute(); //Perform any necessary computations before collision

        void collide(); //Collision step

        void boundaries(); //Boundary calculation

        void initialise(); //Initialisation step

        void computeMomenta(); //Momenta (density, velocity) calculation

        const double& getDensity(int k) const; //Return density at lattice point k

        const std::vector<double>& getVelocity() const; //Return vector of velocity

        const std::vector<double>& getDistribution() const; //Return vector of distribution

    private:

        double computeEquilibrium(const double& orderparam,const double* velocity,
                                  const int idx,const int k) const; //Calculate equilibrium in direction idx with a given
                                                        //density and velocity

        double computeModelForce(int xyz,int k) const; //Calculate forces specific to the model in direction xyz

        double computeForces(int xyz,int k) const; //Calculate other forces in direction xyz

        double computeCollisionQ(double& sum,int k,const double& old,const double& orderparam,
                                 const double* velocity,const int idx) const; //Calculate collision
                                                                                           //at index idx

        double computeOrderParameter(const double* distribution,int k) const; //Calculate the order parameter
                                                                              //corresponding to the relative
                                                                              //concentrations of each phase

        static constexpr double m_Tau=1.0; //TEMPORARY relaxation time

        static constexpr double m_InverseTau=1.0/m_Tau; //TEMPORARY inverse relaxation time

        OrderParameter m_OrderParameter; //Order Parameter

        Velocity m_Velocity; //Velocity

        typename DataType<typename traits::Stencil>::DistributionData& m_Distribution=m_Data.getDistributionObject();
            //Distributions

        DataType<typename traits::Stencil> m_Data; //MOVE THIS TO BASE

        vector<double>& orderparameter=m_OrderParameter.getParameter(); //Reference to vector of order parameters

        vector<double>& velocity=m_Velocity.getParameter(); //Reference to vector of velocities

        vector<double>& distribution=m_Distribution.getDistribution(); //Reference to vector of distributions

        enum{x=0,y=1,z=2}; //Indices corresponding to x, y, z directions

        
        const int& N;
        const int& LX;
        const int& LY;
        const int& LZ;

        double m_Gamma=1;

        ChemicalPotential m_ChemicalPotential;

        GradientOrderParameter m_GradOrderParameter;

        LaplacianOrderParameter m_LaplacianOrderParameter;

        typename traits::Forces mt_Forces; //MOVE THIS TO BASE
        typename traits::Boundaries mt_Boundaries; //MOVE THIS TO BASE

        Geometry m_Geometry; //MOVE THIS TO BASE
        
};

template<int ndim,class traits>
const double& Binary<ndim,traits>::getDensity(int k) const{ //This needs to be renamed

    return orderparameter[k]; //Return reference to order parameter at point k

}

template<int ndim,class traits>
const std::vector<double>& Binary<ndim,traits>::getVelocity() const{

    return velocity; //Return reference to velocity vector

}

template<int ndim,class traits>
const std::vector<double>& Binary<ndim,traits>::getDistribution() const{

    return distribution; //Return reference to distribution vector

}

template<int ndim,class traits>
void Binary<ndim,traits>::precompute(){

    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=LY*LZ*MAXNEIGHBORS;k<N-MAXNEIGHBORS*LY*LZ;k++){ //loop over k

        if constexpr(std::tuple_size<typename traits::Forces::getTupleType>::value!=0){ //Check if there is at least one element
                                                                          //in F
            std::apply([k](auto&... forces){//See Algorithm.hh for explanation of std::apply
                (forces.precompute(k),...);
            }, mt_Forces.getTuple());
        }
        else;
        
    }
    //CHECK DATA TYPE NEEDS THIS
    m_Distribution.getDistribution().swap(m_Distribution.getDistributionOld()); //swap old and new distributions
                                                                                //before collision
}

template<int ndim,class traits>
double Binary<ndim,traits>::computeForces(int xyz,int k) const{

    if constexpr(std::tuple_size<typename traits::Forces::getTupleType>::value!=0){
        return std::apply([xyz,k](auto&... forces){
                return (forces.compute(xyz,k)+...);
            }, mt_Forces.getTuple());
    }
    else return 0;

}

template<int ndim,class traits>
void Binary<ndim,traits>::collide(){

    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=LY*LZ*MAXNEIGHBORS;k<N-MAXNEIGHBORS*LY*LZ;k++){ //loop over k

        double* distribution=m_Distribution.getDistributionPointer(k);
        double* old_distribution=m_Distribution.getDistributionOldPointer(k);
        double sum=0;
        for (int idx=traits::Stencil::Q-1;idx>=0;--idx){ //loop over discrete velocity directions
            //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"
            //std::cout<<m_Distribution.streamIndex(k,idx)<<std::endl;
            m_Distribution.getDistributionPointer(m_Distribution.streamIndex(k,idx))[idx]=computeCollisionQ(sum,k,old_distribution[idx],orderparameter[k],&velocity[k*traits::Stencil::D],idx);
            
        }
        
    }
    #ifdef MPIPARALLEL
    m_Data.communicateDistribution();
    #endif
}

template<int ndim,class traits>
void Binary<ndim,traits>::boundaries(){

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
        else;

    }

    
}

template<int ndim,class traits>
void Binary<ndim,traits>::initialise(){ //Initialise model

    m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    
    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=LY*LZ*MAXNEIGHBORS;k<N-MAXNEIGHBORS*LY*LZ;k++){ //loop over k

        double* distribution=m_Distribution.getDistributionPointer(k);
        double* old_distribution=m_Distribution.getDistributionOldPointer(k);
        m_ChemicalPotential.getParameter(k)=0;
        int xx=computeX(LY,LZ,k);

        if (xx>=LX/2)orderparameter[k]=1.0; //Set order parameter to 1 initially (This will change)
        else orderparameter[k]=-1.0;
        double sum=0;
        for (int idx=traits::Stencil::Q-1;idx>=0;idx--){

            double equilibrium;
            if (idx>0) equilibrium=computeEquilibrium(orderparameter[k],&velocity[k*traits::Stencil::D],idx,k);
            else equilibrium=orderparameter[k]-sum;

            distribution[idx]=equilibrium; //Set distributions to equillibrium
            old_distribution[idx]=equilibrium;        

        }

    }
}


template<int ndim,class traits>
void Binary<ndim,traits>::computeMomenta(){ //Calculate order parameter

    
    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=LY*LZ*MAXNEIGHBORS;k<N-MAXNEIGHBORS*LY*LZ;k++){ //Loop over k
        double* distribution=m_Distribution.getDistributionPointer(k);
        orderparameter[k]=computeOrderParameter(distribution,k);

    }
    #ifdef MPIPARALLEL
    m_Data.communicate(m_OrderParameter);
    #endif
    
}

template<int ndim,class traits>
double Binary<ndim,traits>::computeCollisionQ(double& sum,const int k,const double& old,const double& orderparam,
                                         const double* velocity,const int idx) const{
                                        //Calculate collision step at a given velocity index at point k
    
    double forcexyz[traits::Stencil::D]; //Temporary array storing force in each cartesian direction

    //Force is the sum of model forces and given forces
    for(int xyz=0;xyz< traits::Stencil::D;xyz++) forcexyz[xyz]=computeModelForce(xyz,k)+computeForces(xyz,k);
    
    //Sum of collision + force contributions
    if (idx>0) {
        double eq=CollisionBase<typename traits::Stencil>::collideSRT(old,computeEquilibrium(orderparam,velocity,idx,k),m_InverseTau)
              +CollisionBase<typename traits::Stencil>::forceGuoSRT(forcexyz,velocity,m_InverseTau,idx);
        sum+=eq;
        
        return eq;
    }
    else return m_OrderParameter.getParameter(k)-sum+CollisionBase<typename traits::Stencil>::forceGuoSRT(forcexyz,velocity,m_InverseTau,idx);

}


template<int ndim,class traits>
double Binary<ndim,traits>::computeEquilibrium(const double& orderparam,const double* velocity,const int idx,const int k) const{

    return traits::Stencil::Weights[idx]*(m_ChemicalPotential.getParameter(k)*m_Gamma/traits::Stencil::Cs2+orderparam*CollisionBase<typename traits::Stencil>::computeVelocityFactor(velocity,idx));

}

template<int ndim,class traits>
double Binary<ndim,traits>::computeModelForce(int k,int xyz) const{
    
    return 0;

}


template<int ndim,class traits>
double Binary<ndim,traits>::computeOrderParameter(const double* distribution,int k) const{//Order parameter calculation
    //Order parameter is the sum of distributions plus any source/correction terms
    if constexpr(std::tuple_size<typename traits::Forces::getTupleType>::value!=0){
        return CollisionBase<typename traits::Stencil>::computeZerothMoment(distribution)
        +std::apply([k](auto&... tests){
            return (tests.computeDensitySource(k)+...);
        }, mt_Forces.getTuple());
    }
    //CHANGE THIS SO FIRST/SECOND MOMENT COMPUTATION IS DONE IN DISTRIBUTION
    else return CollisionBase<typename traits::Stencil>::computeZerothMoment(distribution);

}

#endif