#ifndef BINARY_HEADER
#define BINARY_HEADER
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include <utility>

template<class traits>
class Binary:CollisionBase<typename traits::Stencil>{
    public:
        Binary(typename traits::Forces& forces,typename traits::Boundaries& boundaries):mt_Forces(forces),mt_Boundaries(boundaries){
            
        }
        Binary(typename traits::Forces& forces):mt_Forces(forces),mt_Boundaries(std::tuple<>()){
            
        }
        Binary(typename traits::Boundaries& boundaries):mt_Forces(std::tuple<>()),mt_Boundaries(boundaries){
            
        }

        void precompute();

        void collide();

        void boundaries();

        void initialise();

        void computeMomenta();

        const double& getDensity(int k) const;

        const std::vector<double>& getVelocity() const;

        const std::vector<double>& getDistribution() const;

    private:

        double computeEquilibrium(const double& orderparam,const std::vector<double>& velocity, const int idx) const;

        double computeModelForce(int xyz,int k) const;

        double computeForces(int xyz,int k) const;

        double computeCollisionQ(int k,const double& old,const double& orderparam,const std::vector<double>& velocity,const int idx) const;

        double computeOrderParameter(const double* distribution,int k) const;

        static constexpr double m_Tau=1.0;

        static constexpr double m_InverseTau=1.0/m_Tau;

        OrderParameter<double> m_OrderParameter;

        Velocity<double,typename traits::Stencil> m_Velocity;

        Distribution<typename traits::Stencil,Binary<traits>>& m_Distribution=m_Data.template getDistributionObject();

        typename traits::Data<Binary<traits>> m_Data;

        vector<double>& orderparameter=m_OrderParameter.getParameter();

        vector<double>& velocity=m_Velocity.getParameter();

        vector<double>& distribution=m_Distribution.getDistribution();

        enum{x=0,y=1,z=2};

        typename traits::Forces mt_Forces;
        typename traits::Boundaries mt_Boundaries;

        Geometry m_Geometry;
        
};

template<class traits>
const double& Binary<traits>::getDensity(int k) const{

    return orderparameter[k];

}

template<class traits>
const std::vector<double>& Binary<traits>::getVelocity() const{

    return velocity;

}

template<class traits>
const std::vector<double>& Binary<traits>::getDistribution() const{

    return distribution;

}

template<class traits>
void Binary<traits>::precompute(){

    int k=0;

    while(k>=0){

        if constexpr(std::tuple_size<typename traits::Forces>::value!=0){
            std::apply([k](auto&... tests){
                (tests.precompute(k),...);
            }, mt_Forces);
        }
        else;

        k = m_Data.iterateFluid(k);

    }

    m_Distribution.getDistribution().swap(m_Distribution.getDistributionOld());

}

template<class traits>
double Binary<traits>::computeForces(int xyz,int k) const{
    if constexpr (std::tuple_size<typename traits::Forces>::value!=0){
        return std::apply([xyz,k](auto&... tests){
                return (tests.compute(xyz,k)+...);
            }, mt_Forces);
    }
    else return 0;

}

template<class traits>
void Binary<traits>::collide(){

    int k=0;

    while(k>=0){

        double* distribution=m_Distribution.getDistributionPointer(k);
        double* old_distribution=m_Distribution.getDistributionOldPointer(k);

        for (int idx=0;idx<traits::Stencil::Q;idx++){

            m_Distribution.getDistributionPointer(m_Distribution.streamIndex(k,idx))[idx]=computeCollisionQ(k,old_distribution[idx],orderparameter[k],velocity,idx);

        }

        k = m_Data.iterateFluid(k);

    }
}

template<class traits>
void Binary<traits>::boundaries(){

    int k=0;

    while(k>=0){
        if constexpr(std::tuple_size<typename traits::Boundaries>::value!=0){
            std::apply([this,k](auto&... tests){
                (tests.compute(this->m_Distribution,k),...);
            }, mt_Boundaries);
        }
        else;

        k = m_Data.iterateFluid(k);

    }

}

template<class traits>
void Binary<traits>::initialise(){

    m_Data.generateNeighbors();
    
    int k=0;

    while(k>=0){

        double* distribution=m_Distribution.getDistributionPointer(k);
        double* old_distribution=m_Distribution.getDistributionOldPointer(k);

        orderparameter[k]=1.0;

        for (int idx=0;idx<traits::Stencil::Q;idx++){

            double equilibrium=computeEquilibrium(orderparameter[k],velocity,idx);

            distribution[idx]=equilibrium;
            old_distribution[idx]=equilibrium;        

        }

        k = m_Data.iterateFluid(k);

    }
}


template<class traits>
void Binary<traits>::computeMomenta(){

    int k=0;

    while(k>=0){

        double* distribution=m_Distribution.getDistributionPointer(k);

        orderparameter[k]=computeOrderParameter(distribution,k);

        k = m_Data.iterateFluid(k);

    }
}

template<class traits>
double Binary<traits>::computeCollisionQ(const int k,const double& old,const double& orderparam,const std::vector<double>& velocity,const int idx) const{

    std::array<double,traits::Stencil::D> forcexyz;

    for(int xyz=0;xyz< traits::Stencil::D;xyz++) forcexyz[xyz]=computeModelForce(xyz,k)+computeForces(xyz,k);
    
    return CollisionBase<typename traits::Stencil>::collideSRT(old,computeEquilibrium(orderparam,velocity,idx),m_InverseTau)
              +CollisionBase<typename traits::Stencil>::forceSRT(forcexyz,velocity,m_InverseTau,idx);

}


template<class traits>
double Binary<traits>::computeEquilibrium(const double& orderparam,const std::vector<double>& velocity,const int idx) const{

    return orderparam*CollisionBase<typename traits::Stencil>::computeGamma(velocity,idx);

}

template<class traits>
double Binary<traits>::computeModelForce(int k,int xyz) const{

    return 0.0;

}

template<class traits>
double Binary<traits>::computeOrderParameter(const double* distribution,int k) const{
    if constexpr(std::tuple_size<typename traits::Forces>::value!=0){
        return CollisionBase<typename traits::Stencil>::computeFirstMoment(distribution)+std::apply([k](auto&... tests){
            return (tests.computeDensitySource(k)+...);
        }, mt_Forces);
    }
    //CHANGE THIS SO FIRST/SECOND MOMENT COMPUTATION IS DONE IN DISTRIBUTION
    else return CollisionBase<typename traits::Stencil>::computeFirstMoment(distribution);

}

#endif