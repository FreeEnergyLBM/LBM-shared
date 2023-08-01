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
 * There are classes to specify the tau dependence of the collision operator and classes for SRT and MRT collision.
 * The CollisionBase class in this file will be inherited by LBM models to provide basic operations in LBM such as 
 * momenta calculation. If you want to implement a new collision operator, do it here so it can be used by all models.
 */


/**
 * \brief The NoTauDependence class specifies a function to calculate the tau prefactor of the forcing term, which will
 *        just return 1 (as there is no dependence on tau). This can be used as a prefactor in the forcing terms found 
 *        in Forcing.hh.
 */
struct NoTauDependence {

    /**
     * \brief The calculateTauPrefactor function just returns 1.0 in this case, as we want no dependence on tau.
     * \tparam TLattice LatticeProperties class for the system.
     * \param itau 1.0/tau (where tau is the relaxation time).
     * \return 1.0
     */
    template<class TLattice>
    static inline double calculateTauPrefactor(const double& itau) { 

        return 1.;

    } 

};

/**
 * \brief The GuoPrefactor class specifies a function to calculate the tau prefactor of the forcing term, which will 
 *        return the guo forcing prefactor that depends on tau (See [1]). This can be used as a prefactor in the
 *        forcing terms found in Forcing.hh.
 */
struct GuoPrefactor {

    /**
     * \brief The calculateTauPrefactor function will return 1.0-0.5*dt/tau in this case [1].
     * \tparam TLattice LatticeProperties class for the system.
     * \param itau 1.0/tau (where tau is the relaxation time).
     * \return Standard forcing relaxation dependence.
     */
    template<class TLattice>
    static inline double calculateTauPrefactor(const double& itau) {

        return 1.-0.5*TLattice::DT*itau;

    } 

};

/**
 * \brief The SplitGuoPrefactor class specifies a function to calculate the tau prefactor of the forcing term, which will 
 *        return the half of the guo forcing prefactor that depends on tau (See [1]). This can be used as a prefactor in the
 *        forcing terms found in Forcing.hh. This has very niche uses, e.g. for the taehun lee model when you also need non-
 *        equilibrium boundaries.
 */
struct SplitGuoPrefactor {

    /**
     * \brief The calculateTauPrefactor function will return -0.5*dt/tau in this case [1].
     * \tparam TLattice LatticeProperties class for the system.
     * \param itau 1.0/tau (where tau is the relaxation time).
     * \return Split forcing relaxation dependence.
     */
    template<class TLattice>
    static inline double calculateTauPrefactor(const double& itau) { 

        return -0.5*TLattice::DT*itau;

    } 

};

/**
 * \brief The SRT class provides internal functions to calculate the collision operator, Omega, based on the relaxation
 *        time (tau). For SRT (a.k.a. BGK), this is just 1.0/tau. This class also handles the application of forces in the
 *        collision step. Give this to a trait class using SetCollisionOperator<SRT>.
 * \tparam TStencil Velocity stencil (e.g. D2Q9).
 */
template<class TStencil>
class SRT{

    public:

        /**
         * \brief The Omega function will return the SRT/BGK collision operator when given 1.0/tau [1].
         * \param itau 1.0/tau (where tau is the relaxation time).
         * \return 1.0/tau.
         */
        static const inline double& Omega(const double& itau) {return itau;}

        /**
         * \brief The initialise function will perform initialisation for SRT to save time later (currently does nothing) [1].
         * \tparam TLattice LatticeProperties class for the system.
         * \param tau1 Min/Max bound of tau (where tau is the relaxation time).
         * \param tau2 Min/Max bound of tau (alternate to tau1).
         */
        template<class TLattice>
        static inline void initialise(double tau1,double tau2){}; 

        /**
         * \brief The initialise function will perform initialisation for SRT to save time later (currently does nothing) [1].
         * \tparam TLattice LatticeProperties class for the system.
         * \tparam TForceTuple Type of the tuple containing the forces (and source terms) applied to the model.
         * \param forces Tuple containing objects of all the forces (and source terms) applied to the model.
         * \param tau1 Min/Max bound of tau (where tau is the relaxation time).
         * \param tau2 Min/Max bound of tau (alternate to tau1).
         */
        template<class TLattice, typename TForceTuple>
        static inline void initialise(const TForceTuple& forces, double tau1,double tau2){}; 

        /**
         * \brief The collide function returns the post collision distribution [1].
         * \tparam TLattice LatticeProperties class for the system.
         * \param old Pointer to first element of array containing old distributions.
         * \param equilibrium Pointer to first element of array containing updated equilibrium distributions.
         * \param itau Reference to 1.0/tau parameter.
         * \param idx Discrete velocity direction index.
         * \return Updated distribution information.
         */
        template<class TLattice>
        static inline double collide(const double* old, const double* equilibrium, const double& itau, int idx){

            return old[idx] - TLattice::DT * Omega(itau) * (old[idx] - equilibrium[idx]); 

        }

        /**
         * \brief The forcing function will return the forcing contribution to the collision step.
         * \tparam TLattice LatticeProperties class for the system.
         * \tparam TPrefactorType Type of the tau dependant prefactor for the forcing.
         * \tparam TForceTuple Type of the tuple containing the forces (and source terms) applied to the model.
         * \param forces Tuple containing objects of all the forces (and source terms) applied to the model.
         * \param forcearray Pointer to first element of array containing the forces without the tau dependance applied yet.
         * \param itau Reference to 1.0/tau parameter.
         * \param idx Discrete velocity direction index.
         * \return Updated distribution information.
         */
        template<class TLattice, class TPrefactorType, typename TForceTuple>
        static inline double forcing(const TForceTuple& forces, const double* forcearray, const double& itau, int idx){

            return remove_const_and_reference<TPrefactorType>::type::template calculateTauPrefactor<TLattice>(itau) * forcearray[idx]; 

        }

};

/**
 * \brief The MRT class provides internal functions to calculate the collision operator, Omega, based on the relaxation
 *        time (tau). For MRT, each moment has a different relaxation frequency, so this is contained within a matrix.
 *        Information can be found in [1]. This class also handles the application of forces in the collision step using
 *        MRT. Give this to a trait class using SetCollisionOperator<MRT>.
 * \tparam TStencil Velocity stencil (e.g. D2Q9).
 */
template<class TStencil>
class MRT{
    private:

        /**
         * \brief Constuctor is private, so objecs cannot be instantiated outside of this class.
         */
        MRT(){}

        static constexpr int numberofelements = 5000; //!<Number of elements for MRT matrix lookup table
        static constexpr int Qsquared = TStencil::Q * TStencil::Q; //!<Number of velocity directions squared

        double mTaumax; //!<Maximum relaxation time
        double mTaumin; //!<Minimum relaxation time

        double mTauIdxPrefactor; //!<Maximum relaxation time

        double mMRTMatrix[numberofelements * TStencil::Q * TStencil::Q];

        template<typename TForceTuple>
        inline auto& _ForcingMapStore(const TForceTuple& ft){

            static auto ForcingMap = std::apply([this](auto&... forces){

                ct_map<kv<typename remove_const_and_reference<decltype(forces)>::type::Method::Prefactor,
                                std::array<double,this->numberofelements*TStencil::Q*TStencil::Q>>...> tempmap;

                return tempmap;

            }, ft);

            return ForcingMap;

        }

        template<typename TForceTuple>
        inline auto& getForcingMap(const TForceTuple& ft){

            return _ForcingMapStore(ft);

        }

        template<typename TForceTuple>
        inline auto& getForcingMap(const TForceTuple& ft) const{

            return const_cast<typename std::remove_const<decltype(_ForcingMapStore(ft))>::type>(_ForcingMapStore(ft));

        }
        
        template<class TLattice>
        inline void generateMRTTau(double tau, int tauidx, const double (&inverse)[Qsquared]); 

        template<class TLattice, class TTauPrefactor, typename TForceTuple>
        inline void generateMRTTauForcing(const TForceTuple& forces,
                                          double tau, int tauidx, const double (&inverse)[Qsquared]);

        const inline int getTauIdx(double tau) const{

            int tauidx = mTauIdxPrefactor * (tau - mTaumin);
            
            if(tauidx < 0) return 0;
            if(tauidx > numberofelements - 1) return numberofelements - 1;

            return tauidx;

        }

    public:
    
        MRT(MRT const&)             = delete;
        void operator=(MRT const&)  = delete;

        static inline MRT& getInstance(){
            static MRT instance;
            return instance;
        } 

        static const inline double* Omega(const double& itau) {

            static MRT& mrt = getInstance();
            return &(mrt.mMRTMatrix[mrt.getTauIdx(1. / itau) * Qsquared]);

        }

        template<class TTauPrefactor>
        static const inline double* ForcePrefactor(TTauPrefactor& prefactor, const double& itau){

            static MRT& mrt = getInstance();
            return &(prefactor[mrt.getTauIdx(1. / itau) * Qsquared]);

        }

        template<class TLattice>
        static inline void initialise(double tau1, double tau2); 

        template<class TLattice, typename TForceTuple>
        static inline void initialise(const TForceTuple& forces, double tau1,double tau2); 

        template<class TLattice>
        static inline double collide(const double* old, const double* equilibrium, const double& itau, int idx){

            double collisionsum = 0;

            auto MRTArray = Omega(itau);

            for (int sumidx = 0; sumidx<  TStencil::Q; sumidx++){
                collisionsum += TLattice::DT * MRTArray[idx * TStencil::Q + sumidx] 
                                             * (old[sumidx] - equilibrium[sumidx]);
                
            }
            
            return old[idx]-collisionsum;

        }

        template<class TLattice, class prefactortype, typename TForceTuple>
        static inline double forcing(const TForceTuple& forces, const double* forcearray, const double& itau, int idx){

            double forcesum = 0;

            static MRT& mrt = getInstance();
            const auto& prefactor = remove_const_and_reference<decltype(mrt.getForcingMap(forces))>::type::template get<
                                                                    typename std::remove_const<prefactortype>::type>::val;
            auto MRTForcingArray = ForcePrefactor(prefactor, itau);

            for (int sumidx = 0; sumidx < TStencil::Q; sumidx++){

                forcesum += MRTForcingArray[idx * TStencil::Q + sumidx] * (forcearray[sumidx]);

            }

            return forcesum;

        }
         
};

template<class TStencil>
template<class TLattice>
inline void MRT<TStencil>::generateMRTTau(double tau, int tauidx, const double (&Minverse)[MRT<TStencil>::Qsquared]) {

    double weightmatrix[TStencil::Q * TStencil::Q] = {};
    double weightedmoments[TStencil::Q * TStencil::Q] = {};

    for(int i = 0; i < TStencil::Q; i++){

        weightmatrix[i * TStencil::Q + i] = TStencil::MRTWeights(1. / tau)[i];

    }
    
    for(int i = 0; i < TStencil::Q; i++){
        for(int j = 0; j < TStencil::Q; j++){
            for(int ii = 0; ii < TStencil::Q; ii++) {

                weightedmoments[j * TStencil::Q + i] += weightmatrix[j * TStencil::Q + ii] 
                                                        * TStencil::MRTMatrix[ii * TStencil::Q + i];

            }
        }
    }
    
    for(int i = 0; i < TStencil::Q; i++){
        for(int j = 0; j < TStencil::Q; j++){
            for(int ii = 0; ii < TStencil::Q; ii++) { 

                mMRTMatrix[tauidx * TStencil::Q * TStencil::Q 
                                        + j * TStencil::Q + i] += Minverse[j * TStencil::Q + ii]
                                                                  * weightedmoments[ii * TStencil::Q + i];

            }
        }
    }

}

template<class TStencil>
template<class TLattice, class TTauPrefactor, typename TForceTuple>
inline void MRT<TStencil>::generateMRTTauForcing(const TForceTuple& forces,
                                                 double tau, int tauidx,
                                                 const double (&Minverse)[MRT<TStencil>::Qsquared]) {

    if constexpr (remove_const_and_reference< decltype(getInstance().getForcingMap(forces))>::type
                                                        ::template keyexists<TTauPrefactor>::exists){

        double weightmatrixforcing[TStencil::Q * TStencil::Q] = {};
        double weightedmomentsforcing[TStencil::Q * TStencil::Q] = {};

        for(int i = 0; i < TStencil::Q; i++){
        
            weightmatrixforcing[i * TStencil::Q + i] =
                TTauPrefactor::template calculateTauPrefactor<TLattice>(TStencil::MRTWeights(1. / tau)[i]); 

        }
        
        for(int i = 0; i < TStencil::Q; i++){
            for(int j = 0; j < TStencil::Q; j++){
                for(int ii = 0; ii < TStencil::Q; ii++) {  

                    weightedmomentsforcing[j * TStencil::Q + i] +=
                        weightmatrixforcing[j * TStencil::Q + ii]*TStencil::MRTMatrix[ii * TStencil::Q + i];

                }
            }
        }
        using forcingmaptype = typename std::remove_const<
                                    typename std::remove_reference<
                                        decltype(getInstance().getForcingMap(forces))>::type>::type;
        for(int i = 0; i < TStencil::Q; i++){
            for(int j = 0; j < TStencil::Q; j++){

                forcingmaptype::template get<TTauPrefactor>::val[tauidx * TStencil::Q * TStencil::Q 
                                                                              + j * TStencil::Q + i] = 0;
                for(int ii = 0; ii < TStencil::Q; ii++) { 

                    forcingmaptype::template get<TTauPrefactor>
                            ::val[tauidx * TStencil::Q * TStencil::Q
                                               + j * TStencil::Q + i] += Minverse[j * TStencil::Q + ii]
                                                                         * weightedmomentsforcing[ii * TStencil::Q 
                                                                                                               + i];
                    
                }
            }
        }
    }
    
}

template<class TStencil>
template<class TLattice>
inline void MRT<TStencil>::initialise(double tau1, double tau2) {

    #pragma omp master
    {

    getInstance().mTaumax = std::max(tau1, tau2);
    getInstance().mTaumin = std::min(tau1, tau2);
    getInstance().mTauIdxPrefactor = (numberofelements - 1) / (getInstance().mTaumax - getInstance().mTaumin);

    double MomentsInverse[TStencil::Q * TStencil::Q] = {};
    double mag[TStencil::Q] = {};

    for(int j = 0; j < TStencil::Q; j++){ 

        mag[j] = 0; 

        for(int i = 0; i < TStencil::Q; i++){

            mag[j] += TStencil::MRTMatrix[j * TStencil::Q + i]*TStencil::MRTMatrix[j * TStencil::Q + i];

        }
        for(int i = 0; i < TStencil::Q; i++){

            MomentsInverse[i * TStencil::Q + j] = TStencil::MRTMatrix[j * TStencil::Q + i] / mag[j];
            
        }

    }

    for(int tauidx = 0; tauidx < numberofelements; tauidx++){

        double tau = getInstance().mTaumin 
                     + (double)tauidx * (getInstance().mTaumax - getInstance().mTaumin)
                       / ((double)numberofelements - 1.0);

        getInstance().template generateMRTTau<TLattice>(tau, tauidx, MomentsInverse);
        
    }
    
    }

}

template<class TStencil>
template<class TLattice,typename TForceTuple>
inline void MRT<TStencil>::initialise(const TForceTuple& forces, double tau1, double tau2) {

    #pragma omp master
    {

    getInstance().mTaumax = std::max(tau1, tau2);
    getInstance().mTaumin = std::min(tau1, tau2);
    getInstance().mTauIdxPrefactor = (numberofelements - 1) / (getInstance().mTaumax - getInstance().mTaumin);

    double MomentsInverse[TStencil::Q * TStencil::Q] = {};
    double mag[TStencil::Q] = {};

    for(int j = 0; j < TStencil::Q; j++){ 

        mag[j] = 0; 

        for(int i = 0; i < TStencil::Q; i++){

            mag[j] += TStencil::MRTMatrix[j * TStencil::Q + i]*TStencil::MRTMatrix[j * TStencil::Q + i];

        }
        for(int i = 0; i < TStencil::Q; i++){

            MomentsInverse[i * TStencil::Q + j] = TStencil::MRTMatrix[j * TStencil::Q + i] / mag[j];
            
        }

    }

    for(int tauidx = 0; tauidx < numberofelements; tauidx++){

        double tau = getInstance().mTaumin 
                     + (double)tauidx * (getInstance().mTaumax - getInstance().mTaumin)
                       / ((double)numberofelements - 1.0);

        getInstance().template generateMRTTau<TLattice>(tau, tauidx, MomentsInverse);

        std::apply([tau, tauidx, MomentsInverse, forces](auto&... force) {

            (getInstance().template generateMRTTauForcing<TLattice,typename decltype(getMethod(force))::Prefactor>(forces, tau, tauidx, MomentsInverse), ...);

        }, forces);
            
    }
    
    }

}

/**
 * \brief The CollisionBase class provides functions that perform basic LBM calculations e.g. collision operators.
 * This class takes a TStencil as a template argument, as the velocity discretisation information and weights is 
 * needed. The class has public functions for collision terms, momenta calculations, force terms and the common
 * velocity depenence of equilibrium distibutions.
 * \tparam TStencil Velocity TStencil for the model inheriting from this class.
 */
template<class TLattice, class TStencil>
class CollisionBase {
    
    public:
        
        /**
         * \brief computeGamma computes first and second order velocity dependence of the equilibrium distributions, as well as the non velocity dependent part.
         * \param velocity Pointer to velocity vector at the current TLattice point.
         * \param idx The discrete velocity index (e.g. 0-8 for D2Q9).
         * \return 1 + velocity dependence of equilibrium.
         */
        static inline double computeGamma(const double* velocity, int idx);

        /**
         * \brief This will sum the distributions in each direction to calculate the zeroth moment.
         * \param distribution Pointer to distribution vector at the current TLattice point.
         * \return Zeroth moment distributions in each direction.
         */
        static inline double computeZerothMoment(const double *distribution);

        /**
         * \brief This will sum the distributions times the velocity vector
         *        in each direction to calculate the first moment.
         * \param distribution Pointer to distribution vector at the current TLattice point.
         * \param xyz Cartesian direction of to calculate zeroth moment.
         * \return First moment of distributions.
         */
        static inline double computeFirstMoment(const double *distribution, int xyz); 

        /**
         * \brief computeGamma computes first and second order velocity dependence of the equilibrium distributions.
         * \param velocity Pointer to velocity vector at the current TLattice point.
         * \param idx The discrete velocity index (e.g. 0-8 for D2Q9).
         * \return Velocity dependence of equilibrium.
         */
        static inline double computeVelocityFactor(const double* velocity, int idx);
        
    private:

        enum{ x = 0, y = 1, z = 2 };

};

/**
 * \details The computeGamma function will return the standard second order equilibrium distribution divided
 *          by density. This is calcualted as Weights*(1+velocity factor), where "velocity factor" is the velocity
 *          dependence of the equilibrium.
 */
template<class TLattice, class TStencil>
inline double CollisionBase<TLattice,TStencil>::computeGamma(const double* velocity, int idx) {
    
    return TStencil::Weights[idx] * (1.0 + computeVelocityFactor(velocity, idx)); 

};

/**
 * \details This function returns the velocity dependence of the equilibrium distributions. This is seperate
 *          from computeGamma as sometimes this is needed seperately from the usual equilibrium term. First, dot
 *          products of velocity with velocity and velocity with the discrete c_i vectors in the TStencil are
 *          calculated. These are then normalised with respect to the TLattice sound speed and the velocity
 *          factor is returned.
 */
template<class TLattice, class TStencil>
inline double CollisionBase<TLattice,TStencil>::computeVelocityFactor(const double* velocity, int idx) {

    double ci_dot_velocity = (TStencil::Ci_x[idx] * velocity[0]);
    double velocity_dot_velocity = pow(velocity[0], 2);

    if constexpr (TStencil::D > 1) {
        ci_dot_velocity += (TStencil::Ci_y[idx] * velocity[1]); //Dot product of Ci (discrete velocity)
                                                                    //vector and velocity
        velocity_dot_velocity += pow(velocity[1], 2);
    }
    if constexpr (TStencil::D > 2) {
        ci_dot_velocity += (TStencil::Ci_z[idx] * velocity[2]); //Dot product of Ci (discrete velocity)
                                                                    //vector and velocity
        velocity_dot_velocity += pow(velocity[2], 2);
    }

    return (ci_dot_velocity) / TStencil::Cs2
           + (ci_dot_velocity * ci_dot_velocity) / (2.0 * TStencil::Cs2 * TStencil::Cs2) //Return velocity factor
           - (velocity_dot_velocity) / (2.0 * TStencil::Cs2);

};

/**
 * \details This function returns the zeroth moment of the distributions. This is just the sum of distributions
 *          in each discrete direction (so the sum over 9 directions for D2Q9);
 */
template<class TLattice, class TStencil>
inline double CollisionBase<TLattice,TStencil>::computeZerothMoment(const double *distribution) {

    double zerothmoment = 0;

    for (int idx = 0; idx <TStencil::Q; idx++) {

        zerothmoment += distribution[idx]; //Sum distribution over Q

    }

    return zerothmoment; //And return the sum

}

/**
 * \details This function returns the first moment of the distributions. This is the sum of the distributions
 *          multiplied by the TStencil velocity vector c_i for each i in the choesn cartesian direction.
 */
template<class TLattice, class TStencil>
inline double CollisionBase<TLattice,TStencil>::computeFirstMoment(const double *distribution, int xyz) {

    double firstmoment = 0;

    for (int idx = 0; idx <TStencil::Q; idx++){
        
        firstmoment += (distribution[idx]*TStencil::Ci_xyz(xyz)[idx]); //Sum distribution times Ci over Q
        
    }
    
    return firstmoment; //Return first moment corresponding to velocity in given direction ("xyz")

}