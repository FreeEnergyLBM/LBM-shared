#pragma once
#include <string>
#include <memory>
#include <utility>
#include <vector>
#include <iostream>
#include "Lattice.hh"
#include "Global.hh"
#include "Stencil.hh"

/**
 * \file Collide.hh
 * \brief Contains base class with commonly used functions for the collision and momentum calculation steps in LBM.
 * The class in this file will be inherited by LBM models to provide basic operations in LBM such as SRT collision
 * or Guo forcing. If you want to implement a new collision or forcing operator, do it here so it can be used by
 * all models.
 */

/**
 * \brief The CollisionBase class provides functions that perform basic LBM calculations e.g. collision operators.
 * This class takes a stencil as a template argument, as the velocity discretisation information and weights is 
 * needed. The class has public functions for collision terms, momenta calculations, force terms and the common
 * velocity depenence of equilibrium distibutions.
 * \tparam stencil Velocity stencil for the model inheriting from this class.
 */
template<class lattice, class stencil>
class CollisionBase {

    static_assert(std::is_base_of<Stencil, stencil>(), "ERROR: invalid stencil specified in traits class.");
    static_assert(stencil::D == lattice::m_NDIM, "ERROR: The chosen stencil must match the number of lattice dimensions in the lattice properties.");
    
    public:
        
        /**
         * \brief computeGamma computes first and second order velocity dependence of the equilibrium distributions, as well as the non velocity dependent part.
         * \param velocity Pointer to velocity vector at the current lattice point.
         * \param idx The discrete velocity index (e.g. 0-8 for D2Q9).
         * \return 1 + velocity dependence of equilibrium.
         */
        inline double computeGamma(const double* velocity, const int idx) const;

        /**
         * \brief This will sum the distributions in each direction to calculate the zeroth moment.
         * \param distribution Pointer to distribution vector at the current lattice point.
         * \return Zeroth moment distributions in each direction.
         */
        inline double computeZerothMoment(const double *distribution) const;

        /**
         * \brief This will sum the distributions times the velocity vector
         *        in each direction to calculate the first moment.
         * \param distribution Pointer to distribution vector at the current lattice point.
         * \param xyz Cartesian direction of to calculate zeroth moment.
         * \return First moment of distributions.
         */
        inline double computeFirstMoment(const double *distribution, const int xyz) const; 

        /**
         * \brief This will compute the single relaxation time (BGK) collision step.
         * \param old Distribution at previous timestep.
         * \param equilibrium Equilibrium distribution at the current timestep.
         * \param tau Relaxation time (chosen to provide a desired viscosity).
         * \return Post collision distribution.
         */
        inline double collideSRT(const double& old, const double& equilibrium, const double& tau) const;

        /**
         * \brief This will compute the forcing term using Guo forcing.
         * \param force Array containing the total force in the cartesian directions.
         * \param velocity Pointer to velocity vector at the current lattice point.
         * \param itau Inverse relaxation time (1.0/tau) (chosen to provide a desired viscosity).
         * \param idx The discrete velocity index (e.g. 0-8 for D2Q9).
         * \return Forcing term in chosen velocity index direction.
         */
        inline double forceGuoSRT(const double force[stencil::D ],const double* velocity,
                            const double& itau, const int idx) const;
        
        /**
         * \brief computeGamma computes first and second order velocity dependence of the equilibrium distributions.
         * \param velocity Pointer to velocity vector at the current lattice point.
         * \param idx The discrete velocity index (e.g. 0-8 for D2Q9).
         * \return Velocity dependence of equilibrium.
         */
        inline double computeVelocityFactor(const double* velocity, const int idx) const;
        
    private:

        enum{ x = 0, y = 1, z = 2 };
        
        //static constexpr auto& ma_Weights = stencil::Weights;

        static constexpr double m_Cs2 = stencil::Cs2;

};

/**
 * \details The computeGamma function will return the standard second order equilibrium distribution divided
 *          by density. This is calcualted as Weights*(1+velocity factor), where "velocity factor" is the velocity
 *          dependence of the equilibrium.
 */
template<class lattice, class stencil>
inline double CollisionBase<lattice,stencil>::computeGamma(const double* velocity, const int idx) const {
    
    return stencil::Weights[idx] * (1.0 + computeVelocityFactor(velocity, idx)); 

};

/**
 * \details This function returns the velocity dependence of the equilibrium distributions. This is seperate
 *          from computeGamma as sometimes this is needed seperately from the usual equilibrium term. First, dot
 *          products of velocity with velocity and velocity with the discrete c_i vectors in the stencil are
 *          calculated. These are then normalised with respect to the lattice sound speed and the velocity
 *          factor is returned.
 */
template<class lattice, class stencil>
inline double CollisionBase<lattice,stencil>::computeVelocityFactor(const double* velocity, const int idx) const {
    
    double ci_dot_velocity = 0;
    double velocity_dot_velocity = 0;

    for (int xyz = 0; xyz <stencil::D; xyz++) {

        ci_dot_velocity += (stencil::Ci_xyz(xyz)[idx] * velocity[xyz]); //Dot product of Ci (discrete velocity)
                                                                    //vector and velocity
        velocity_dot_velocity += (velocity[xyz] * velocity[xyz]); //Dot product of velocity and velocity

    }

    return (ci_dot_velocity) / m_Cs2
           + (ci_dot_velocity * ci_dot_velocity) / (2.0 * m_Cs2 * m_Cs2) //Return velocity factor
           - (velocity_dot_velocity) / (2.0 * m_Cs2);

};

/**
 * \details This function returns the zeroth moment of the distributions. This is just the sum of distributions
 *          in each discrete direction (so the sum over 9 directions for D2Q9);
 */
template<class lattice, class stencil>
inline double CollisionBase<lattice,stencil>::computeZerothMoment(const double *distribution) const {

    double zerothmoment = 0;

    for (int idx = 0; idx <stencil::Q; idx++) {

        zerothmoment += distribution[idx]; //Sum distribution over Q

    }

    return zerothmoment; //And return the sum

}

/**
 * \details This function returns the first moment of the distributions. This is the sum of the distributions
 *          multiplied by the stencil velocity vector c_i for each i in the choesn cartesian direction.
 */
template<class lattice, class stencil>
inline double CollisionBase<lattice,stencil>::computeFirstMoment(const double *distribution,const int xyz) const {

    double firstmoment = 0;

    for (int idx = 0; idx <stencil::Q; idx++){

        firstmoment+=(distribution[idx]*stencil::Ci_xyz(xyz)[idx]); //Sum distribution times Ci over Q
        
    }
    
    return firstmoment; //Return first moment corresponding to velocity in given direction ("xyz")

}

/**
 * \details This computes the SRT/BGK collision step. This is simply the old distribution in each direction minus
 *          the difference between the old and equilibrium distributions divided by the relaxation time, tau.
 */
template<class lattice, class stencil>
inline double CollisionBase<lattice,stencil>::collideSRT(const double& old, const double& equilibrium, const double& itau) const{

    return old - lattice::m_DT * itau * (old - equilibrium); //SRT colision step. Old corresponds to the old distribution and itau is
                                       //the inverse of the relaxation time

}

/**
 * \details This computes the Guo forcing for the SRT collision operator. The dot product of velocity with the
 *          stencil velocity vectors is calculated and then used to calculate the forcing term.
 */
template<class lattice, class stencil>
inline double CollisionBase<lattice,stencil>::forceGuoSRT(const double force[stencil::D],
                                        const double* velocity, const double& itau,
                                        const int idx) const { //Guo forcing

    double ci_dot_velocity = 0;
    double forceterm = 0;
    double prefactor = (1 - lattice::m_DT * itau / 2.0) * stencil::Weights[idx]; //Prefactor for Guo forcing

    for (int xyz=0;xyz<stencil::D;xyz++){

        ci_dot_velocity += (stencil::Ci_xyz(xyz)[idx] * velocity[xyz]); //Dot product of discrete velocity vector
                                                                    //with velocity
    }
    for (int xyz = 0; xyz <stencil::D; xyz++) {

        forceterm += prefactor * (((stencil::Ci_xyz(xyz)[idx] - velocity[xyz]) / m_Cs2
                               + ci_dot_velocity * stencil::Ci_xyz(xyz)[idx] / (m_Cs2 * m_Cs2)) * force[xyz]); //Force
                                                                                                     //Calculation
    }
    
    return forceterm;

}
