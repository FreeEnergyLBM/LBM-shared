#ifndef BINARY_HEADER
#define BINARY_HEADER
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include <utility>

//Binary.hh: Contains the details of the LBM model to solve an equation for phase separation. Each
//Model is given a "traits" class that contains stencil, data, force and boundary information

template<class traits>
class Binary:CollisionBase<typename traits::Stencil>{ //Inherit from base class to avoid repetition of common
                                                      //calculations
    public:
        //Constructors to construct tuples of forces and boundaries
        Binary(typename traits::Forces& forces,typename traits::Boundaries& boundaries):mt_Forces(forces),mt_Boundaries(boundaries){
            
        }
        Binary(typename traits::Forces& forces):mt_Forces(forces),mt_Boundaries(std::tuple<>()){
            
        }
        Binary(typename traits::Boundaries& boundaries):mt_Forces(std::tuple<>()),mt_Boundaries(boundaries){
            
        }
        Binary():mt_Forces(std::tuple<>()),mt_Boundaries(std::tuple<>()){
            
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
                                  const int idx) const; //Calculate equilibrium in direction idx with a given
                                                        //density and velocity

        double computeModelForce(int xyz,int k) const; //Calculate forces specific to the model in direction xyz

        double computeForces(int xyz,int k) const; //Calculate other forces in direction xyz

        double computeCollisionQ(int k,const double& old,const double& orderparam,
                                 const double* velocity,const int idx) const; //Calculate collision
                                                                                           //at index idx

        double computeOrderParameter(const double* distribution,int k) const; //Calculate the order parameter
                                                                              //corresponding to the relative
                                                                              //concentrations of each phase

        static constexpr double m_Tau=1.0; //TEMPORARY relaxation time

        static constexpr double m_InverseTau=1.0/m_Tau; //TEMPORARY inverse relaxation time

        OrderParameter<double> m_OrderParameter; //Order Parameter

        Velocity<double,typename traits::Stencil> m_Velocity; //Velocity

        Distribution_Base<typename traits::Stencil>& m_Distribution=m_Data.template getDistributionObject();
            //Distributions

        typename traits::Data m_Data; //MOVE THIS TO BASE

        vector<double>& orderparameter=m_OrderParameter.getParameter(); //Reference to vector of order parameters

        vector<double>& velocity=m_Velocity.getParameter(); //Reference to vector of velocities

        vector<double>& distribution=m_Distribution.getDistribution(); //Reference to vector of distributions

        enum{x=0,y=1,z=2}; //Indices corresponding to x, y, z directions

        typename traits::Forces mt_Forces; //MOVE THIS TO BASE
        typename traits::Boundaries mt_Boundaries; //MOVE THIS TO BASE

        Geometry m_Geometry; //MOVE THIS TO BASE
        
};

template<class traits>
const double& Binary<traits>::getDensity(int k) const{ //This needs to be renamed

    return orderparameter[k]; //Return reference to order parameter at point k

}

template<class traits>
const std::vector<double>& Binary<traits>::getVelocity() const{

    return velocity; //Return reference to velocity vector

}

template<class traits>
const std::vector<double>& Binary<traits>::getDistribution() const{

    return distribution; //Return reference to distribution vector

}

template<class traits>
void Binary<traits>::precompute(){

    int k=LY*LZ*MAXNEIGHBORS;
    k = m_Data.iterateFluid0(k,false);
    while(k>=0){ //loop over k

        if constexpr(std::tuple_size<typename traits::Forces>::value!=0){ //Check if there is at least one element
                                                                          //in F
            std::apply([k](auto&... forces){//See Algorithm.hh for explanation of std::apply
                (forces.precompute(k),...);
            }, mt_Forces);
        }
        else;

        k = m_Data.iterateFluid(k,false); //increment k
        
    }

    m_Distribution.getDistribution().swap(m_Distribution.getDistributionOld()); //swap old and new distributions
                                                                                //before collision
}

template<class traits>
double Binary<traits>::computeForces(int xyz,int k) const{

    if constexpr(std::tuple_size<typename traits::Forces>::value!=0){
        return std::apply([xyz,k](auto&... forces){
                return (forces.compute(xyz,k)+...);
            }, mt_Forces);
    }
    else return 0;

}

template<class traits>
void Binary<traits>::collide(){

    int k=LY*LZ*MAXNEIGHBORS;
    k = m_Data.iterateFluid0(k,false);
    while(k>=0){ //loop over k

        double* distribution=m_Distribution.getDistributionPointer(k);
        double* old_distribution=m_Distribution.getDistributionOldPointer(k);

        for (int idx=0;idx<traits::Stencil::Q;idx++){ //loop over discrete velocity directions
            //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"
            m_Distribution.getDistributionPointer(m_Distribution.streamIndex(k,idx))[idx]=computeCollisionQ(k,old_distribution[idx],orderparameter[k],&velocity[k*traits::Stencil::D],idx);
            
        }
        
        k = m_Data.iterateFluid(k,false); //increment k

        
    }

    m_Data.communicateDistribution();

}

template<class traits>
void Binary<traits>::boundaries(){

    int k=0;
    k = m_Data.iterateSolid0(k,true);
    while(k>=0){ //loop over k

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
        
        k = m_Data.iterateSolid(k,true); //increment k

    }
    
}

template<class traits>
void Binary<traits>::initialise(){ //Initialise model

    m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    
    int k=LY*LZ*MAXNEIGHBORS;
    k = m_Data.iterateFluid0(k,false);
    while(k>=0){ //loop over k

        double* distribution=m_Distribution.getDistributionPointer(k);
        double* old_distribution=m_Distribution.getDistributionOldPointer(k);

        orderparameter[k]=1.0; //Set order parameter to 1 initially (This will change)

        for (int idx=0;idx<traits::Stencil::Q;idx++){

            double equilibrium=computeEquilibrium(orderparameter[k],&velocity[k*traits::Stencil::D],idx);

            distribution[idx]=equilibrium; //Set distributions to equillibrium
            old_distribution[idx]=equilibrium;        

        }

        k = m_Data.iterateFluid(k,false); //increment k

    }
}


template<class traits>
void Binary<traits>::computeMomenta(){ //Calculate order parameter

    int k=LY*LZ*MAXNEIGHBORS;
    k = m_Data.iterateFluid0(k,false);
    while(k>=0){ //Loop over k

        double* distribution=m_Distribution.getDistributionPointer(k);

        orderparameter[k]=computeOrderParameter(distribution,k);

        k = m_Data.iterateFluid(k,false); //increment k

    }

    m_Data.communicate(m_OrderParameter);

}

template<class traits>
double Binary<traits>::computeCollisionQ(const int k,const double& old,const double& orderparam,
                                         const double* velocity,const int idx) const{
                                        //Calculate collision step at a given velocity index at point k

    std::array<double,traits::Stencil::D> forcexyz; //Temporary array storing force in each cartesian direction

    //Force is the sum of model forces and given forces
    for(int xyz=0;xyz< traits::Stencil::D;xyz++) forcexyz[xyz]=computeModelForce(xyz,k)+computeForces(xyz,k);
    
    //Sum of collision + force contributions
    return CollisionBase<typename traits::Stencil>::collideSRT(old,computeEquilibrium(orderparam,velocity,idx),m_InverseTau)
              +CollisionBase<typename traits::Stencil>::forceSRT(forcexyz,velocity,m_InverseTau,idx);

}


template<class traits>
double Binary<traits>::computeEquilibrium(const double& orderparam,const double* velocity,const int idx) const{

    return orderparam*CollisionBase<typename traits::Stencil>::computeGamma(velocity,idx);//Equilibrium is density
                                                                                          //times gamma in this
                                                                                          //case (THIS IS WRONG)

}

template<class traits>
double Binary<traits>::computeModelForce(int k,int xyz) const{

    return 0.0; //No model force in this case

}

template<class traits>
double Binary<traits>::computeOrderParameter(const double* distribution,int k) const{//Order parameter calculation
    //Order parameter is the sum of distributions plus any source/correction terms
    if constexpr(std::tuple_size<typename traits::Forces>::value!=0){
        return CollisionBase<typename traits::Stencil>::computeFirstMoment(distribution)+std::apply([k](auto&... tests){
            return (tests.computeDensitySource(k)+...);
        }, mt_Forces);
    }
    //CHANGE THIS SO FIRST/SECOND MOMENT COMPUTATION IS DONE IN DISTRIBUTION
    else return CollisionBase<typename traits::Stencil>::computeFirstMoment(distribution);

}

#endif