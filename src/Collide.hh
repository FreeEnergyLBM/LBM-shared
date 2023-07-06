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

template<class stencil>
class SRT{
    public:
        static const inline double& Omega(const double& itau) {return itau;}
        template<class lattice>
        static const inline double ForcePrefactor(const double& itau) {return 1.-0.5*lattice::DT*itau;}
        template<class lattice>
        static inline void initialise(double tau1,double tau2){}; 
        template<class lattice>
        static inline double collide(const double* old, const double* equilibrium, const double& itau, int idx){
            return old[idx] - lattice::DT * Omega(itau) * (old[idx] - equilibrium[idx]); 
        }
        template<class lattice>
        static inline double forcing(const double* forces, const double& itau, int idx){
            return ForcePrefactor<lattice>(itau)*forces[idx]; 
        }
};

template<class stencil>
class MRT{
    private:
        static constexpr int numberofelements = 5000;
        double m_Taumax;
        double m_Taumin;
        
        double m_MRTMatrix[numberofelements*stencil::Q*stencil::Q];
        double m_MRTMatrixForcing[numberofelements*stencil::Q*stencil::Q];
        template<class lattice>
        inline void generateMRTTau(double tau,int tauidx,const double (&inverse)[stencil::Q*stencil::Q]); 
        MRT(){}
        inline int getTauIdx(double tau){
            int tauidx = (numberofelements-1)*(tau-m_Taumin)
            /(m_Taumax-m_Taumin);
            
            if(tauidx < 0) tauidx = 0; 
            if(tauidx > numberofelements-1) tauidx = numberofelements-1;
            //if(tau>0.6&&tau<0.9)std::cout<<tau<<" "<<tauidx<<std::endl;
            return tauidx;
        }
    public:
        static inline MRT& getInstance(){
            static MRT instance;
            return instance;
        } 
        MRT(MRT const&)             = delete;
        void operator=(MRT const&)  = delete;
        static const inline double* Omega(const double& itau) {
            return &(getInstance().m_MRTMatrix[getInstance().getTauIdx(1./itau)*stencil::Q*stencil::Q]);
        }
        static const inline double* ForcePrefactor(const double& itau){
            return &(getInstance().m_MRTMatrixForcing[getInstance().getTauIdx(1./itau)*stencil::Q*stencil::Q]);
        }
        template<class lattice>
        static inline void initialise(double tau1,double tau2);  

        template<class lattice>
        static inline double collide(const double* old, const double* equilibrium, const double& itau, int idx){
            double collisionsum=0;
            for (int sumidx=0; sumidx<stencil::Q; sumidx++){
                collisionsum+=lattice::DT * getInstance().Omega(itau)[idx*stencil::Q+sumidx] * (old[sumidx] - equilibrium[sumidx]);
                
            }
            
            return old[idx]-collisionsum;
        }

        template<class lattice>
        static inline double forcing(const double* forces, const double& itau, int idx){
            double forcesum=0;
            for (int sumidx=0; sumidx<stencil::Q; sumidx++){
                forcesum+=getInstance().template ForcePrefactor(itau)[idx*stencil::Q+sumidx] * (forces[sumidx]);
            }
            return forcesum;
        }
         
};

template<class stencil>
template<class lattice>
inline void MRT<stencil>::generateMRTTau(double tau,int tauidx,const double (&Minverse)[stencil::Q*stencil::Q]) {

    double weightmatrix[stencil::Q*stencil::Q] = {};
    double weightmatrixforcing[stencil::Q*stencil::Q] = {};
    double weightedmoments[stencil::Q*stencil::Q] = {};
    double weightedmomentsforcing[stencil::Q*stencil::Q] = {};

    for(int i = 0; i < stencil::Q; i++){
        weightmatrix[i*stencil::Q+i] = stencil::MRTWeights(1./tau)[i];
        weightmatrixforcing[i*stencil::Q+i] = 1.0-lattice::DT*0.5*weightmatrix[i*stencil::Q+i]; 
    }
    
    for(int i = 0; i < stencil::Q; i++){
        for(int j = 0; j < stencil::Q; j++){
            
            for(int ii = 0; ii < stencil::Q; ii++) {
                weightedmoments[j*stencil::Q+i] += weightmatrix[j*stencil::Q+ii]*stencil::MRTMatrix[ii*stencil::Q+i];
                
                weightedmomentsforcing[j*stencil::Q+i] += weightmatrixforcing[j*stencil::Q+ii]*stencil::MRTMatrix[ii*stencil::Q+i];
            }
            
            
        }
    }
    
    for(int i = 0; i < stencil::Q; i++){
        for(int j = 0; j < stencil::Q; j++){
            for(int ii = 0; ii < stencil::Q; ii++) { 
                m_MRTMatrix[tauidx*stencil::Q*stencil::Q+j*stencil::Q+i] += Minverse[j*stencil::Q+ii]*weightedmoments[ii*stencil::Q+i];
                m_MRTMatrixForcing[tauidx*stencil::Q*stencil::Q+j*stencil::Q+i] += Minverse[j*stencil::Q+ii]*weightedmomentsforcing[ii*stencil::Q+i];
            }
            
        }
    }

}

template<class stencil>
template<class lattice>
inline void MRT<stencil>::initialise(double tau1,double tau2) {
    #pragma omp master
    {

    getInstance().m_Taumax=std::max(tau1,tau2);
    getInstance().m_Taumin=std::min(tau1,tau2);

    double MomentsInverse[stencil::Q*stencil::Q] = {};
    double mag[stencil::Q] = {};

    for(int j = 0; j < stencil::Q; j++){ 
        mag[j] = 0; 
        for(int i = 0; i < stencil::Q; i++){
            mag[j] += stencil::MRTMatrix[j*stencil::Q+i]*stencil::MRTMatrix[j*stencil::Q+i];
        }
        for(int i = 0; i < stencil::Q; i++){
            MomentsInverse[i*stencil::Q+j] = stencil::MRTMatrix[j*stencil::Q+i]/mag[j];
            //std::cout<<stencil::MRTMatrix[j*stencil::Q+i]<<std::endl;
            
        }
    }

    for(int tauidx = 0; tauidx < numberofelements; tauidx++){

        double tau = getInstance().m_Taumin + (double)tauidx*(getInstance().m_Taumax-getInstance().m_Taumin)/((double)numberofelements-1.0);

        getInstance().template generateMRTTau<lattice>(tau,tauidx,MomentsInverse);
        

    }
    
    }

}

/**
 * \brief The CollisionBase class provides functions that perform basic LBM calculations e.g. collision operators.
 * This class takes a stencil as a template argument, as the velocity discretisation information and weights is 
 * needed. The class has public functions for collision terms, momenta calculations, force terms and the common
 * velocity depenence of equilibrium distibutions.
 * \tparam stencil Velocity stencil for the model inheriting from this class.
 */
template<class lattice, class stencil>
class CollisionBase {
    
    public:
        
        /**
         * \brief computeGamma computes first and second order velocity dependence of the equilibrium distributions, as well as the non velocity dependent part.
         * \param velocity Pointer to velocity vector at the current lattice point.
         * \param idx The discrete velocity index (e.g. 0-8 for D2Q9).
         * \return 1 + velocity dependence of equilibrium.
         */
        static inline double computeGamma(const double* velocity, const int idx);

        /**
         * \brief This will sum the distributions in each direction to calculate the zeroth moment.
         * \param distribution Pointer to distribution vector at the current lattice point.
         * \return Zeroth moment distributions in each direction.
         */
        static inline double computeZerothMoment(const double *distribution);

        /**
         * \brief This will sum the distributions times the velocity vector
         *        in each direction to calculate the first moment.
         * \param distribution Pointer to distribution vector at the current lattice point.
         * \param xyz Cartesian direction of to calculate zeroth moment.
         * \return First moment of distributions.
         */
        static inline double computeFirstMoment(const double *distribution, const int xyz); 

        /**
         * \brief computeGamma computes first and second order velocity dependence of the equilibrium distributions.
         * \param velocity Pointer to velocity vector at the current lattice point.
         * \param idx The discrete velocity index (e.g. 0-8 for D2Q9).
         * \return Velocity dependence of equilibrium.
         */
        static inline double computeVelocityFactor(const double* velocity, const int idx);
        
    private:

        enum{ x = 0, y = 1, z = 2 };
        
        //static constexpr auto& ma_Weights = stencil::Weights;

};

/**
 * \details The computeGamma function will return the standard second order equilibrium distribution divided
 *          by density. This is calcualted as Weights*(1+velocity factor), where "velocity factor" is the velocity
 *          dependence of the equilibrium.
 */
template<class lattice, class stencil>
inline double CollisionBase<lattice,stencil>::computeGamma(const double* velocity, const int idx) {
    
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
inline double CollisionBase<lattice,stencil>::computeVelocityFactor(const double* velocity, const int idx) {
    
    double ci_dot_velocity = 0;
    double velocity_dot_velocity = 0;

    for (int xyz = 0; xyz <stencil::D; xyz++) {

        ci_dot_velocity += (stencil::Ci_xyz(xyz)[idx] * velocity[xyz]); //Dot product of Ci (discrete velocity)
                                                                    //vector and velocity
        velocity_dot_velocity += (velocity[xyz] * velocity[xyz]); //Dot product of velocity and velocity

    }

    return (ci_dot_velocity) / stencil::Cs2
           + (ci_dot_velocity * ci_dot_velocity) / (2.0 * stencil::Cs2 * stencil::Cs2) //Return velocity factor
           - (velocity_dot_velocity) / (2.0 * stencil::Cs2);

};

/**
 * \details This function returns the zeroth moment of the distributions. This is just the sum of distributions
 *          in each discrete direction (so the sum over 9 directions for D2Q9);
 */
template<class lattice, class stencil>
inline double CollisionBase<lattice,stencil>::computeZerothMoment(const double *distribution) {

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
inline double CollisionBase<lattice,stencil>::computeFirstMoment(const double *distribution,const int xyz) {

    double firstmoment = 0;

    for (int idx = 0; idx <stencil::Q; idx++){

        firstmoment+=(distribution[idx]*stencil::Ci_xyz(xyz)[idx]); //Sum distribution times Ci over Q
        
    }
    
    return firstmoment; //Return first moment corresponding to velocity in given direction ("xyz")

}