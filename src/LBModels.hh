#include "Collide.hh"

template<class stencil,class... forces>
class SingleComponent:CollisionBase<stencil>{
    public:
        void precompute() const;
        double collide(double old) const;
        void initialise() const;
    private:
        double computeEquilibrium() const;
        double computeModelForce() const;

        double computeDensity(const double* distribution) const;
        double computeVelocity() const;
        const double m_Tau=1.0;
        const double m_InverseTau=1.0/m_Tau;
        enum{x=0,y=1,z=2};
        Density m_Density;
        Velocity m_Velocity;
        Distribution1 m_Distribution;
};



template<class stencil,class... forces>
void SingleComponent<stencil,forces...>::precompute() const{
    (forces::precompute(),...);
}

template<class stencil,class... forces>
template<int velocityindex,int ...xyz>
double SingleComponent<stencil,forces...>::collide(double old) const{
    return old+CollisionBase<stencil>::collideSRT(old,computeEquilibrium<velocityindex,xyz...>(*m_Density.getParameter(),m_Velocity.getParameter()),m_InverseTau)+CollisionBase<stencil>::forceSRT<velocityindex,xyz...>(computeModelForce(),velocity,m_InverseTau)+(CollisionBase<stencil>::forceSRT<velocityindex,xyz...>(forces::compute(),velocity,m_InverseTau)+...);
}

template<class stencil,class... forces>
void SingleComponent<stencil,forces...>::initialise() const{
    
}


template<class stencil,class... forces>
template<int i,int ...xyz>
void SingleComponent<stencil,forces...>::computeEquilibrium(const double& density,const double *velocity) const{
    return density*CollisionBase<stencil>::computeGamma<i,xyz...>(velocity);
}

template<class stencil,class... forces>
void SingleComponent<stencil,forces...>::computeModelForce(const double& density,const double *velocity) const{
    return 0.0;
}

template<class stencil,class... forces>
void SingleComponent<stencil,forces...>::computeDensity(const double* distribution) const{
    return computeFirstMoment(distribution,make_index_sequence<stencil::m_Q>)+((forces::computeDensitySource())+...);
}

template<class stencil,class... forces>
template<int d>
void SingleComponent<stencil,forces...>::computeVelocity(const double& density) const{
    return computeSecondMoment<d>(distribution,make_index_sequence<stencil::m_Q>)+((forces::computeVelocitySource())+...);
}