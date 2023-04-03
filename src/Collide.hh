#ifndef COLLIDE_HEADER
#define COLLIDE_HEADER
#include <string>
#include <memory>
#include <utility>
#include <vector>
#include <iostream>
#include "Global.hh"

//Collide.hh: Contains base class with commonly used functions for the collision and momentum calculation
//steps in LBM

template<class stencil> //Stencil is needed for moment calculations
class CollisionBase{
    static_assert(stencil::D==NDIM,"ERROR: The chosen stencil must match the number of lattice dimensions (NDIM) chosen in Global.hh.");
    public:

        double computeGamma(const double* velocity, const int idx) const; //Gamma is the standard
                                                                                       //equilibrium calculation
                                                                                       //divided by density

        double computeFirstMoment(const double *distribution) const; //Sum distributions over Q to calculate
                                                                     //first moment

        double computeSecondMoment(const double *distribution, const int xyz) const; //Sum distributions*C_i
                                                                                     //over Q to calculate
                                                                                     //second moment

        double collideSRT(const double& old,const double& equilibrium,const double& tau) const; //SRT collision
                                                                                                //step

        double forceSRT(const std::array<double,stencil::D> force,const double* velocity,
                        const double& itau,const int idx) const; //SRT force calculation

        double computeVelocityFactor(const double* velocity, const int idx) const; //Second
                                                                                       //order velocity dependence
                                                                                       //of the equilibrium
                                                                                       //distributions times

    private:

        enum{x=0,y=1,z=2}; //Indices of x, y, z directions

        static constexpr auto& ma_Weights=stencil::Weights; //Reference to stencil weights to shorten code
                                                            //somewhat

        static constexpr double m_Cs2=stencil::Cs2; //Again, just to shorten code
        
};

template<class stencil>
double CollisionBase<stencil>::computeGamma(const double* velocity, const int idx) const{

    return ma_Weights[idx]*(1.0+computeVelocityFactor(velocity,idx)); //Weights*(1+velocity factor)
                                                                      //This is the standard equilibrium in
                                                                      //LBM but not multiplied by density
                                                                      //SEE LITERATURE

};

template<class stencil>
double CollisionBase<stencil>::computeVelocityFactor(const double* velocity, const int idx) const{
    //Sometimes the velocity part of the equilibrium is needed seperately so we do this here
    double ci_dot_velocity=0;
    double velocity_dot_velocity=0;

    for (int xyz=0;xyz<stencil::D;xyz++){
        ci_dot_velocity+=(stencil::Ci_xyz(xyz)[idx]*velocity[xyz]); //Dot product of Ci (discrete velocity)
                                                                    //vector and velocity
        velocity_dot_velocity+=(velocity[xyz]*velocity[xyz]); //Dot product of velocity and velocity
    }

    return (ci_dot_velocity)/m_Cs2
           +(ci_dot_velocity*ci_dot_velocity)/(2.0*m_Cs2*m_Cs2) //Return velocity factor
           -(velocity_dot_velocity)/(2.0*m_Cs2);

};

template<class stencil>
double CollisionBase<stencil>::computeFirstMoment(const double *distribution) const{

    double firstmoment=0;

    for (int idx=0;idx<stencil::Q;idx++){
        firstmoment+=distribution[idx]; //Sum distribution over Q
    }

    return firstmoment; //And return

}

template<class stencil>

double CollisionBase<stencil>::computeSecondMoment(const double *distribution,const int xyz) const{

    double secondmoment=0;

    for (int idx=0;idx<stencil::Q;idx++){
        secondmoment+=(distribution[idx]*stencil::Ci_xyz(xyz)[idx]); //Sum distribution times Ci over Q
    }
    
    return secondmoment; //Return second moment corresponding to velocity in given direction ("xyz")

}

template<class stencil>
double CollisionBase<stencil>::collideSRT(const double& old,const double& equilibrium,const double& itau) const{

    return old-itau*(old-equilibrium); //SRT colision step. Old corresponds to the old distribution and itau is
                                       //the inverse of the relaxation time

}

template<class stencil>
double CollisionBase<stencil>::forceSRT(const std::array<double,stencil::D> force,
                                        const double* velocity,const double& itau,
                                        const int idx) const{ //Guo forcing //SEE LITERATURE

    double ci_dot_velocity=0;
    double forceterm=0;
    double prefactor=(1-DT*itau/2.0)*ma_Weights[idx]; //Prefactor for Guo forcing

    for (int xyz=0;xyz<stencil::D;xyz++){
        ci_dot_velocity+=(stencil::Ci_xyz(xyz)[idx]*velocity[xyz]); //Dot product of discrete velocity vector
                                                                    //with velocity
    }
    for (int xyz=0;xyz<stencil::D;xyz++){
        forceterm+=prefactor*(((stencil::Ci_xyz(xyz)[idx]-velocity[xyz])/m_Cs2
                               +ci_dot_velocity*stencil::Ci_xyz(xyz)[idx]/(m_Cs2*m_Cs2))*force[xyz]); //Force
                                                                                                     //Calculation
    }
    
    return forceterm;

}
#endif