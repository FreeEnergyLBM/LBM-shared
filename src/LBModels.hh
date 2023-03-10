#ifndef LBMODELS_HEADER
#define LBMODELS_HEADER
#include "Collide.hh"
#include "Parameters.hh"
#include "Data.hh"
#include <utility>

template<class stencil,template<typename givenstencil> class data,class... forces>
class SingleComponent:CollisionBase<stencil>{
    public:
        SingleComponent(forces&... Forces):mt_Forces(Forces...){
            
        }

        void precompute();

        void collide();

        void initialise();

        void computeMomenta();

        const double& getDensity(int k) const;
        const std::vector<double>& getVelocity() const;
        const std::vector<double>& getDistribution() const;

    private:

        double computeEquilibrium(const double& density,const std::vector<double>& velocity, const int idx) const;

        double computeModelForce(int xyz,int k) const;

        double computeForces(int xyz,int k) const;

        double computeCollisionQ(int k,const double& old,const double& density,const std::vector<double>& velocity,const int idx) const;

        double computeDensity(const double* distribution,int k) const;

        double computeVelocity(const double* distribution,const double& density, const int xyz,int k) const;

        static constexpr double m_Tau=0.75;

        static constexpr double m_InverseTau=1.0/m_Tau;

        Density<double> m_Density;

        Velocity<double,stencil> m_Velocity;

        Distribution<stencil,1> m_Distribution=m_Data.template getDistributionObject<1>();

        vector<double>& density=m_Density.getParameter();
        vector<double>& velocity=m_Velocity.getParameter();
        vector<double>& distribution=m_Distribution.getDistribution();
        enum{x=0,y=1,z=2};
        std::tuple<forces&...> mt_Forces;
        data<stencil> m_Data;
};

template<class stencil,template<typename givenstencil> class data,class... forces>
const double& SingleComponent<stencil,data,forces...>::getDensity(int k) const{
    return density[k];
}

template<class stencil,template<typename givenstencil> class data,class... forces>
const std::vector<double>& SingleComponent<stencil,data,forces...>::getVelocity() const{
    return velocity;
}

template<class stencil,template<typename givenstencil> class data,class... forces>
const std::vector<double>& SingleComponent<stencil,data,forces...>::getDistribution() const{
    return distribution;
}

template<class stencil,template<typename givenstencil> class data,class... forces>
void SingleComponent<stencil,data,forces...>::precompute(){
    int k=0;
    while(k>=0){
        if constexpr (sizeof...(forces)!=0) (std::get<forces&>(mt_Forces).precompute(k),...);
        else;
        k = m_Data.iterate(k);

    }
    m_Distribution.getDistribution().swap(m_Distribution.getDistributionOld());
}

template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeForces(int xyz,int k) const{
    if constexpr (sizeof...(forces)!=0) return (std::get<forces&>(mt_Forces).compute(xyz,k)+...);
    else return 0;
}

template<class stencil,template<typename givenstencil> class data,class... forces>
void SingleComponent<stencil,data,forces...>::collide(){
    int k=0;
    while(k>=0){
        double* distribution=m_Distribution.getDistributionPointer(k);
        double* old_distribution=m_Distribution.getDistributionOldPointer(k);

        for (int idx=0;idx<stencil::Q;idx++){
            distribution[m_Distribution.streamIndex(k,idx)]=computeCollisionQ(k,old_distribution[idx],density[k],velocity,idx);
        }
        k = m_Data.iterate(k);

    }
}

template<class stencil,template<typename givenstencil> class data,class... forces>
void SingleComponent<stencil,data,forces...>::initialise(){
    int k=0;
    while(k>=0){
        double* distribution=m_Distribution.getDistributionPointer(k);
        double* old_distribution=m_Distribution.getDistributionOldPointer(k);
        density[k]=1.0;
        velocity[k*stencil::D+x]=0.0;
        velocity[k*stencil::D+y]=0;
        velocity[k*stencil::D+z]=0;
        for (int idx=0;idx<stencil::Q;idx++){
            double equilibrium=computeEquilibrium(density[k],velocity,idx);
            //std::cout<<computeEquilibrium(density[k],velocity,idx)<<std::endl;

            distribution[idx]=equilibrium;
            //std::cout<<distribution[idx]<<std::endl;
            old_distribution[idx]=equilibrium;        
        }
        k = m_Data.iterate(k);

    }
}


template<class stencil,template<typename givenstencil> class data,class... forces>
void SingleComponent<stencil,data,forces...>::computeMomenta(){
    int k=0;
    while(k>=0){

        double* distribution=m_Distribution.getDistributionPointer(k);
        density[k]=computeDensity(distribution,k);
        velocity[k*stencil::D+x]=computeVelocity(distribution,density[k],x,k);
        velocity[k*stencil::D+y]=computeVelocity(distribution,density[k],y,k);
        velocity[k*stencil::D+z]=computeVelocity(distribution,density[k],z,k);
        k = m_Data.iterate(k);

    }
}

template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeCollisionQ(const int k,const double& old,const double& density,const std::vector<double>& velocity,const int idx) const{
    std::array<double,stencil::D> forcexyz;
    for(int xyz=0;xyz<stencil::D;xyz++) forcexyz[xyz]=computeModelForce(xyz,k)+computeForces(xyz,k);
    
    return CollisionBase<stencil>::collideSRT(old,computeEquilibrium(density,velocity,idx),m_InverseTau)
              +CollisionBase<stencil>::forceSRT(forcexyz,velocity,m_InverseTau,idx);

}


template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeEquilibrium(const double& density,const std::vector<double>& velocity,const int idx) const{
    return density*CollisionBase<stencil>::computeGamma(velocity,idx);
}

template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeModelForce(int k,int xyz) const{
    return 0.0;
}

template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeDensity(const double* distribution,int k) const{
    if constexpr (sizeof...(forces)!=0) return CollisionBase<stencil>::computeFirstMoment(distribution)+((std::get<forces&>(mt_Forces).computeDensitySource(k))+...); //CHANGE THIS SO FIRST/SECOND MOMENT COMPUTATION IS DONE IN DISTRIBUTION
    else return CollisionBase<stencil>::computeFirstMoment(distribution);
}

template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeVelocity(const double* distribution,const double& density, const int xyz,int k) const{
    if constexpr (sizeof...(forces)!=0) return CollisionBase<stencil>::computeSecondMoment(distribution,xyz)+((std::get<forces&>(mt_Forces).computeVelocitySource(xyz,k))+...);
    else return CollisionBase<stencil>::computeSecondMoment(distribution,xyz);
}
#endif