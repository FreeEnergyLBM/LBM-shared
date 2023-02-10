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

        const double& getDensity() const;
        const std::vector<double>& getVelocity() const;
        const std::vector<double>& getDistribution() const;

    private:

        double computeEquilibrium(const double& density,const std::vector<double>& velocity, const int idx) const;

        double computeModelForce(int xyz) const;

        double computeForces(int xyz) const;

        double computeCollisionQ(const double& old,const double& density,const std::vector<double>& velocity,const int idx) const;

        double computeDensity(const double* distribution) const;

        double computeVelocity(const double* distribution,const double& density, const int xyz) const;

        static constexpr double m_Tau=1.0;

        static constexpr double m_InverseTau=1.0/m_Tau;

        Density<double> m_Density;

        Velocity<double,stencil> m_Velocity;

        Distribution<stencil,data,1> m_Distribution;

        vector<double>& density=m_Density.getParameter();
        vector<double>& velocity=m_Velocity.getParameter();
        vector<double>& distribution=m_Distribution.getDistribution();
        enum{x=0,y=1,z=2};
        std::tuple<forces&...> mt_Forces;
};

template<class stencil,template<typename givenstencil> class data,class... forces>
const double& SingleComponent<stencil,data,forces...>::getDensity() const{
    return density[0];
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
    if constexpr (sizeof...(forces)!=0) (std::get<forces&>(mt_Forces).precompute(),...);
    else;
    std::swap(m_Distribution.getDistribution(0),m_Distribution.getDistribution(0));
}

template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeForces(int xyz) const{
    if constexpr (sizeof...(forces)!=0) return (std::get<forces&>(mt_Forces).compute(xyz)+...);
    else return 0;
}

template<class stencil,template<typename givenstencil> class data,class... forces>
void SingleComponent<stencil,data,forces...>::collide(){

    double* distribution=m_Distribution.getDistributionPointer(0);
    double* old_distribution=m_Distribution.getDistributionOldPointer(0);

    for (int idx=0;idx<stencil::Q;idx++){
        distribution[idx]=computeCollisionQ(old_distribution[idx],density[0],velocity,idx);
    }

}

template<class stencil,template<typename givenstencil> class data,class... forces>
void SingleComponent<stencil,data,forces...>::initialise(){

    double* distribution=m_Distribution.getDistributionPointer(0);
    double* old_distribution=m_Distribution.getDistributionOldPointer(0);

    for (int idx=0;idx<stencil::Q;idx++){
        double equilibrium=computeEquilibrium(density[0],velocity,idx);
        distribution[idx]=equilibrium;
        old_distribution[idx]=equilibrium;
        
    }
    density[0]=1.0;
    velocity[x]=0.01;
    velocity[y]=0;
    velocity[z]=0;

}


template<class stencil,template<typename givenstencil> class data,class... forces>
void SingleComponent<stencil,data,forces...>::computeMomenta(){

    double* distribution=m_Distribution.getDistributionPointer(0);
    density[0]=computeDensity(distribution);
    velocity[x]=computeVelocity(distribution,density[0],x);
    velocity[y]=computeVelocity(distribution,density[0],y);
    velocity[z]=computeVelocity(distribution,density[0],z);

}

template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeCollisionQ(const double& old,const double& density,const std::vector<double>& velocity,const int idx) const{
    std::array<double,stencil::D> forcexyz;
    for(int xyz=0;xyz<stencil::D;xyz++) forcexyz[xyz]=computeModelForce(xyz)+computeForces(xyz);
    
    return old+CollisionBase<stencil>::collideSRT(old,computeEquilibrium(density,velocity,idx),m_InverseTau)
              +CollisionBase<stencil>::forceSRT(forcexyz,velocity,m_InverseTau,idx);

}


template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeEquilibrium(const double& density,const std::vector<double>& velocity,const int idx) const{
    return density*CollisionBase<stencil>::computeGamma(velocity,idx);
}

template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeModelForce(int xyz) const{
    return 0.0;
}

template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeDensity(const double* distribution) const{
    if constexpr (sizeof...(forces)!=0) CollisionBase<stencil>::computeFirstMoment(distribution)+((std::get<forces&>(mt_Forces).computeDensitySource())+...); //CHANGE THIS SO FIRST/SECOND MOMENT COMPUTATION IS DONE IN DISTRIBUTION
    else return CollisionBase<stencil>::computeFirstMoment(distribution);
}

template<class stencil,template<typename givenstencil> class data,class... forces>
double SingleComponent<stencil,data,forces...>::computeVelocity(const double* distribution,const double& density, const int xyz) const{
    if constexpr (sizeof...(forces)!=0) return CollisionBase<stencil>::computeSecondMoment(distribution,xyz)+((std::get<forces&>(mt_Forces).computeVelocitySource(xyz))+...);
    else return CollisionBase<stencil>::computeSecondMoment(distribution,xyz);
}
#endif