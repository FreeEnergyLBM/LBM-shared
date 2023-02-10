#ifndef COLLIDE_HEADER
#define COLLIDE_HEADER
#include <string>
#include <memory>
#include <utility>
#include <vector>
#include <iostream>
#include "Global.hh"


template<class stencil>
class CollisionBase{
    public:

        double computeGamma(const std::vector<double>& velocity, const int idx) const;

        double computeFirstMoment(const double *distribution) const;

        double computeSecondMoment(const double *distribution, const int xyz) const;

        double collideSRT(const double& old,const double& equilibrium,const double& tau) const;

        double forceSRT(const std::array<double,stencil::D> force,const std::vector<double>& velocity,const double& itau,const int idx) const;

    private:

        double computeVelocityFactor(const std::vector<double>& velocity, const int idx) const;

        enum{x=0,y=1,z=2};

        static constexpr auto& ma_Weights=stencil::Weights;

        static constexpr double m_Cs2=stencil::Cs2;
        
};

template<class stencil>
double CollisionBase<stencil>::computeGamma(const std::vector<double>& velocity, const int idx) const{

    return ma_Weights[idx]*(1.0+computeVelocityFactor(velocity,idx));

};

template<class stencil>
double CollisionBase<stencil>::computeVelocityFactor(const std::vector<double>& velocity, const int idx) const{

    double ci_dot_velocity=0;
    double velocity_dot_velocity=0;

    for (int xyz=0;xyz<stencil::D;xyz++){
        ci_dot_velocity+=(stencil::Ci_xyz(xyz)[idx]*velocity[xyz]);
        velocity_dot_velocity+=(velocity[xyz]*velocity[xyz]);
    }

    return (ci_dot_velocity)/m_Cs2
           +(ci_dot_velocity*ci_dot_velocity)/(2.0*m_Cs2*m_Cs2)
           -(velocity_dot_velocity)/(2.0*m_Cs2);

};

template<class stencil>
double CollisionBase<stencil>::computeFirstMoment(const double *distribution) const{

    double firstmoment=0;

    for (int idx=0;idx<stencil::Q;idx++){
        firstmoment+=distribution[idx];
    }

    return firstmoment;

}

template<class stencil>

double CollisionBase<stencil>::computeSecondMoment(const double *distribution,const int xyz) const{

    double secondmoment=0;

    for (int idx=0;idx<stencil::Q;idx++){
        secondmoment+=(distribution[idx]*stencil::Ci_xyz(xyz)[idx]);
    }

    return secondmoment;

}

template<class stencil>
double CollisionBase<stencil>::collideSRT(const double& old,const double& equilibrium,const double& itau) const{

    return old-itau*(old-equilibrium);

}

template<class stencil>
double CollisionBase<stencil>::forceSRT(const std::array<double,stencil::D> force,const std::vector<double>& velocity,const double& itau,const int idx) const{

    double ci_dot_velocity=0;
    double forceterm=0;
    double prefactor=(1-DT*itau/2.0)*ma_Weights[idx];

    for (int xyz=0;xyz<stencil::D;xyz++){
        ci_dot_velocity+=(stencil::Ci_xyz(xyz)[idx]*velocity[xyz]);
    }
    for (int xyz=0;xyz<stencil::D;xyz++){
        forceterm+=prefactor*(((stencil::Ci_xyz(xyz)[idx]-velocity[xyz])/m_Cs2+ci_dot_velocity*stencil::Ci_xyz(xyz)[idx]/(m_Cs2*m_Cs2))*force[xyz]);
    }
    
    return forceterm;

}
#endif