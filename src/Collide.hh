#pragma once
#include <string>
#include <memory>
#include <utility>
#include <vector>
#include <iostream>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
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

struct SplitGuoPrefactor {
    template<class TLattice>
    static inline double calculateTauPrefactor(const double& itau) { 
        return -0.5*TLattice::DT*itau;
    } 
};

template<class TStencil>
class SRT{
    public:
        static const inline double& Omega(const double& itau) {return itau;}
        template<class TLattice>
        static inline void initialise(double tau1,double tau2){}; 
        template<class TLattice, typename TForceTuple>
        static inline void initialise(const TForceTuple& forces, double tau1,double tau2){}; 
        template<class TLattice>
        static inline double collide(const double* old, const double* equilibrium, const double& itau, int idx){
            return old[idx] - TLattice::DT * Omega(itau) * (old[idx] - equilibrium[idx]); 
        }
        template<class TLattice, class TPrefactorType, typename TForceTuple>
        static inline double forcing(const TForceTuple& forces, const double* forcearray, const double& itau, int idx){
            //std::cout<<forces[idx]<<std::endl;
            return std::remove_const<typename std::remove_reference<TPrefactorType>::type>::type::template calculateTauPrefactor<TLattice>(itau)*forcearray[idx]; 
        }
};

template<class TStencil>
class MRT{
    private:
        static constexpr int numberofelements = 5000;
        double mTaumax;
        double mTaumin;
        double mTauIdxPrefactor;
        static constexpr int Qsquared = TStencil::Q*TStencil::Q;
        double mMRTMatrix[numberofelements*TStencil::Q*TStencil::Q];
        template<int i>
        void test(){}
        template<typename TForceTuple>
        inline auto& _ForcingMapStore(const TForceTuple& ft){

            static auto ForcingMap = std::apply([this](auto&... forces){//See Algorithm.hh for explanation of std::apply

                ct_map<kv<typename std::remove_const<decltype(std::remove_const<typename std::remove_reference<decltype(forces)>::type>::type::Method::Prefactor)>::type,std::array<double,this->numberofelements*TStencil::Q*TStencil::Q>>...> tempmap;
                //test<ct_map<kv<typename std::remove_const<decltype(std::remove_const<typename std::remove_reference<decltype(forces)>::type>::type::Method::Prefactor)>::type,std::array<double,this->numberofelements*TStencil::Q*TStencil::Q>>...>>();
                return tempmap;

            }, ft);
            return ForcingMap;//ForcingMap;

        }
        template<typename TForceTuple>
        inline auto& getForcingMap(const TForceTuple& ft){

            //static auto& ForcingMap = _ForcingMapStore(ft);
            return _ForcingMapStore(ft);

        }
        template<typename TForceTuple>
        inline auto& getForcingMap(const TForceTuple& ft) const{

            //static auto ForcingMap = std::apply([this](auto&... forces){//See Algorithm.hh for explanation of std::apply

            //    ct_map<kv<typename std::remove_const<decltype(std::remove_reference<decltype(forces)>::type::Method::Prefactor)>::type,std::array<double,this->numberofelements*Qsquared>>...> tempmap;

            //    return tempmap;

            //}, ft);
            //return ForcingMap;

            //static auto& ForcingMap = _ForcingMapStore(ft);
            return const_cast<typename std::remove_const<decltype(_ForcingMapStore(ft))>::type>(_ForcingMapStore(ft));//ForcingMap;

        }
        //std::unordered_map<std::type_index, std::array<double,numberofelements*TStencil::Q*TStencil::Q>> mM_MRTMatrixForcing;
        
        template<class TLattice>
        inline void generateMRTTau(double tau,int tauidx,const double (&inverse)[Qsquared]); 
        template<class TLattice,class TTauPrefactor,typename TForceTuple>
        inline void generateMRTTauForcing(const TTauPrefactor& ForcePrefactor, const TForceTuple& forces, double tau,int tauidx,const double (&inverse)[Qsquared]); 
        MRT(){}
        const inline int getTauIdx(double tau) const{

            int tauidx = mTauIdxPrefactor*(tau-mTaumin);
            
            if(tauidx < 0) return 0;
            if(tauidx > numberofelements-1) return numberofelements-1;
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
            static MRT& mrt=getInstance();
            return &(mrt.mMRTMatrix[mrt.getTauIdx(1./itau)*Qsquared]);
        }
        template<class TTauPrefactor>
        static const inline double* ForcePrefactor(TTauPrefactor& prefactor, const double& itau){
            //std::cout<<typeid(TTauPrefactor).name()<<std::endl;
            static MRT& mrt=getInstance();
            
            return &(prefactor[mrt.getTauIdx(1./itau)*Qsquared]);
        }
        template<class TLattice>
        static inline void initialise(double tau1,double tau2);  
        template<class TLattice,typename TForceTuple>
        static inline void initialise(const TForceTuple& forces, double tau1,double tau2); 

        template<class TLattice>
        static inline double collide(const double* old, const double* equilibrium, const double& itau, int idx){
            double collisionsum=0;
            auto MRTArray = Omega(itau);
            for (int sumidx=0; sumidx<TStencil::Q; sumidx++){
                collisionsum+=TLattice::DT * MRTArray[idx*TStencil::Q+sumidx] * (old[sumidx] - equilibrium[sumidx]);
                
            }
            
            return old[idx]-collisionsum;
        }

        template<class TLattice, class prefactortype,typename TForceTuple>
        static inline double forcing(const TForceTuple& forces, const double* forcearray, const double& itau, int idx){
            double forcesum=0;
            static MRT& mrt=getInstance();
            const auto& prefactor = std::remove_const<typename std::remove_reference<decltype(mrt.getForcingMap(forces))>::type>::type::template get<typename std::remove_const<prefactortype>::type>::val;
            auto MRTForcingArray = ForcePrefactor(prefactor,itau);
            for (int sumidx=0; sumidx<TStencil::Q; sumidx++){
                forcesum+=MRTForcingArray[idx*TStencil::Q+sumidx] * (forcearray[sumidx]);
                //std::cout<<forces[sumidx]<<std::endl;
            }
            return forcesum;
        }
         
};

template<class TStencil>
template<class TLattice>
inline void MRT<TStencil>::generateMRTTau(double tau,int tauidx,const double (&Minverse)[MRT<TStencil>::Qsquared]) {

    double weightmatrix[TStencil::Q*TStencil::Q] = {};
    double weightedmoments[TStencil::Q*TStencil::Q] = {};

    for(int i = 0; i < TStencil::Q; i++){
        weightmatrix[i*TStencil::Q+i] = TStencil::MRTWeights(1./tau)[i];
    }
    
    for(int i = 0; i < TStencil::Q; i++){
        for(int j = 0; j < TStencil::Q; j++){
            
            for(int ii = 0; ii < TStencil::Q; ii++) {
                weightedmoments[j*TStencil::Q+i] += weightmatrix[j*TStencil::Q+ii]*TStencil::MRTMatrix[ii*TStencil::Q+i];
            }
            
        }
    }
    
    for(int i = 0; i < TStencil::Q; i++){
        for(int j = 0; j < TStencil::Q; j++){
            for(int ii = 0; ii < TStencil::Q; ii++) { 
                mMRTMatrix[tauidx*TStencil::Q*TStencil::Q+j*TStencil::Q+i] += Minverse[j*TStencil::Q+ii]*weightedmoments[ii*TStencil::Q+i];
            }
        }
    }

}

template<class TStencil>
template<class TLattice, class TTauPrefactor, typename TForceTuple>
inline void MRT<TStencil>::generateMRTTauForcing(const TTauPrefactor& ForcePrefactor, const TForceTuple& forces, double tau,int tauidx,const double (&Minverse)[MRT<TStencil>::Qsquared]) {
    //std::cout<<"init "<<typeid(TTauPrefactor).name()<<std::endl;
    if constexpr (std::remove_const<typename std::remove_reference<decltype(getInstance().getForcingMap(forces))>::type>::type::template keyexists<TTauPrefactor>::exists){
        double weightmatrixforcing[TStencil::Q*TStencil::Q] = {};
        double weightedmomentsforcing[TStencil::Q*TStencil::Q] = {};

        for(int i = 0; i < TStencil::Q; i++){
            weightmatrixforcing[i*TStencil::Q+i] = ForcePrefactor.template calculateTauPrefactor<TLattice>(TStencil::MRTWeights(1./tau)[i]); 
        }
        
        for(int i = 0; i < TStencil::Q; i++){
            for(int j = 0; j < TStencil::Q; j++){
                
                for(int ii = 0; ii < TStencil::Q; ii++) {                
                    weightedmomentsforcing[j*TStencil::Q+i] += weightmatrixforcing[j*TStencil::Q+ii]*TStencil::MRTMatrix[ii*TStencil::Q+i];
                }
                
            }
        }
        
        for(int i = 0; i < TStencil::Q; i++){
            for(int j = 0; j < TStencil::Q; j++){
                std::remove_const<typename std::remove_reference<decltype(getInstance().getForcingMap(forces))>::type>::type::template get<TTauPrefactor>::val[tauidx*TStencil::Q*TStencil::Q+j*TStencil::Q+i] = 0;
                for(int ii = 0; ii < TStencil::Q; ii++) { 
                    std::remove_const<typename std::remove_reference<decltype(getInstance().getForcingMap(forces))>::type>::type::template get<TTauPrefactor>::val[tauidx*TStencil::Q*TStencil::Q+j*TStencil::Q+i] += Minverse[j*TStencil::Q+ii]*weightedmomentsforcing[ii*TStencil::Q+i];
                    
                }
            }
        }
    }
    
}

template<class TStencil>
template<class TLattice>
inline void MRT<TStencil>::initialise(double tau1,double tau2) {
    #pragma omp master
    {

    getInstance().mTaumax=std::max(tau1,tau2);
    getInstance().mTaumin=std::min(tau1,tau2);
    getInstance().mTauIdxPrefactor=(numberofelements-1)/(getInstance().mTaumax-getInstance().mTaumin);
    double MomentsInverse[TStencil::Q*TStencil::Q] = {};
    double mag[TStencil::Q] = {};

    for(int j = 0; j < TStencil::Q; j++){ 
        mag[j] = 0; 
        for(int i = 0; i < TStencil::Q; i++){
            mag[j] += TStencil::MRTMatrix[j*TStencil::Q+i]*TStencil::MRTMatrix[j*TStencil::Q+i];
        }
        for(int i = 0; i < TStencil::Q; i++){
            MomentsInverse[i*TStencil::Q+j] = TStencil::MRTMatrix[j*TStencil::Q+i]/mag[j];
            //std::cout<<TStencil::MRTMatrix[j*TStencil::Q+i]<<std::endl;
            
        }
    }

    for(int tauidx = 0; tauidx < numberofelements; tauidx++){

        double tau = getInstance().mTaumin + (double)tauidx*(getInstance().mTaumax-getInstance().mTaumin)/((double)numberofelements-1.0);

        getInstance().template generateMRTTau<TLattice>(tau,tauidx,MomentsInverse);
        

    }
    
    }

}

template<class TStencil>
template<class TLattice,typename TForceTuple>
inline void MRT<TStencil>::initialise(const TForceTuple& forces, double tau1,double tau2) {
    #pragma omp master
    {

    getInstance().mTaumax=std::max(tau1,tau2);
    getInstance().mTaumin=std::min(tau1,tau2);
    getInstance().mTauIdxPrefactor=(numberofelements-1)/(getInstance().mTaumax-getInstance().mTaumin);
    double MomentsInverse[TStencil::Q*TStencil::Q] = {};
    double mag[TStencil::Q] = {};

    for(int j = 0; j < TStencil::Q; j++){ 
        mag[j] = 0; 
        for(int i = 0; i < TStencil::Q; i++){
            mag[j] += TStencil::MRTMatrix[j*TStencil::Q+i]*TStencil::MRTMatrix[j*TStencil::Q+i];
        }
        for(int i = 0; i < TStencil::Q; i++){
            MomentsInverse[i*TStencil::Q+j] = TStencil::MRTMatrix[j*TStencil::Q+i]/mag[j];
            //std::cout<<TStencil::MRTMatrix[j*TStencil::Q+i]<<std::endl;
            
        }
    }

    for(int tauidx = 0; tauidx < numberofelements; tauidx++){

        double tau = getInstance().mTaumin + (double)tauidx*(getInstance().mTaumax-getInstance().mTaumin)/((double)numberofelements-1.0);

        getInstance().template generateMRTTau<TLattice>(tau,tauidx,MomentsInverse);
        std::apply([tau,tauidx,MomentsInverse,forces](auto&... force) {

                (getInstance().template generateMRTTauForcing<TLattice>(decltype(getMethod(force))::Prefactor,forces,tau,tauidx,MomentsInverse),...);

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
        static inline double computeGamma(const double* velocity, const int idx);

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
        static inline double computeFirstMoment(const double *distribution, const int xyz); 

        /**
         * \brief computeGamma computes first and second order velocity dependence of the equilibrium distributions.
         * \param velocity Pointer to velocity vector at the current TLattice point.
         * \param idx The discrete velocity index (e.g. 0-8 for D2Q9).
         * \return Velocity dependence of equilibrium.
         */
        static inline double computeVelocityFactor(const double* velocity, const int idx);
        
    private:

        enum{ x = 0, y = 1, z = 2 };
        
        //static constexpr auto& ma_Weights = TStencil::Weights;

};

/**
 * \details The computeGamma function will return the standard second order equilibrium distribution divided
 *          by density. This is calcualted as Weights*(1+velocity factor), where "velocity factor" is the velocity
 *          dependence of the equilibrium.
 */
template<class TLattice, class TStencil>
inline double CollisionBase<TLattice,TStencil>::computeGamma(const double* velocity, const int idx) {
    
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
inline double CollisionBase<TLattice,TStencil>::computeVelocityFactor(const double* velocity, const int idx) {

    double ci_dot_velocity = (TStencil::Ci_x[idx] * velocity[0]);
    double velocity_dot_velocity = pow(velocity[0],2);

    if constexpr (TStencil::D>1) {
        ci_dot_velocity += (TStencil::Ci_y[idx] * velocity[1]); //Dot product of Ci (discrete velocity)
                                                                    //vector and velocity
        velocity_dot_velocity += pow(velocity[1],2);
    }
    if constexpr (TStencil::D>2) {
        ci_dot_velocity += (TStencil::Ci_z[idx] * velocity[2]); //Dot product of Ci (discrete velocity)
                                                                    //vector and velocity
        velocity_dot_velocity += pow(velocity[2],2);
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
inline double CollisionBase<TLattice,TStencil>::computeFirstMoment(const double *distribution,const int xyz) {

    double firstmoment = 0;

    for (int idx = 0; idx <TStencil::Q; idx++){
        
        firstmoment+=(distribution[idx]*TStencil::Ci_xyz(xyz)[idx]); //Sum distribution times Ci over Q
        
    }
    
    return firstmoment; //Return first moment corresponding to velocity in given direction ("xyz")

}