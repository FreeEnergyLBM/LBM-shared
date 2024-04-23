#pragma once
#include <vector>
#include<string>
#include<map>
#include<any>
#include "Global.hh"
#include <fstream>
#include <iostream>

//Parameters.hh: This file details how macroscopic quantities are stored and interacted with.
//The Parameter class contains a static vector and functions to return the parameters stored.
//The reason for this static vector is to ensure that, if I use "Density" in multiple
//classes, the value will be consistent between classes. Note that, for every template configuration, a new
//static variable is created, so I pass each class to itself as a template to ensure the parameters are unique.


template<class TObj, class TLattice, typename T, int TNum=1> //obj template will guarantee a unique instance of the class with its own
                                       //static vector. I pass the class to itself to guarantee this
                                       //
                                       //T determines the type stored in the vector, TNum
                                       //determines the number of directions for each lattice point (eg you might
                                       //want two directions corresponding to velocity in x and velocity in y at
                                       //each point.
class Parameter {

    public:

        std::map<int,bool> mmInitialised; //Change this to std::set

        static constexpr int mNum=TNum;
        static constexpr int mInstances=TObj::instances;
        static constexpr int mDirections=mNum/TObj::instances;
        static constexpr const char *mName = TObj::mName;

        using ParamType = T;

        inline void save(int t);

        std::vector<T> mv_Parameter; //Static vector (Does not change between objects of the class)
        std::vector<T> mv_CommParameter; //Parameter vector reordered for convenience of communication

        template<class, class, int>
        friend class ParameterSingleton;

        inline std::vector<T>& getParameter() { //Get a vector containing the parameters in the communication region
            return mv_Parameter;
        }

        inline std::vector<T>& getCommParameter() { //Get a vector containing the parameters in the communication region
            return mv_CommParameter;
        }

    private:

        Parameter() {

            mv_Parameter.resize(TNum * TLattice::N); //Resize to the desired size
            mv_CommParameter.resize(TNum * 4 * TLattice::HaloSize); //Resize to the size of the communication region

        }
        
        Parameter(const Parameter& other) {}
   
};

template<class TObj, typename T=double, int TNumPrefactor=1>
class ParameterSingleton {

    public:
        static constexpr int instances=TNumPrefactor;
        using ParamType = T;

        template<int num=1>
        constexpr static auto idxlambda = [](int idx){return (instances>1) ? idx : (num>1)*idx;};
        
        template<class TLattice, int TNum=1>
        static inline Parameter<TObj,TLattice,T,TNumPrefactor*TNum>& getInstance() { static Parameter<TObj,TLattice,T,TNumPrefactor*TNum> instance;
                                                                          return instance;};

        template<class TLattice, int TNum=1>
        static inline std::vector<T>& get() { //Returns vector containing the parameter

            //static Parameter<TObj,TLattice,T,TNumPrefactor*TNum>& instance = getInstance<TLattice,TNum>();
            return getInstance<TLattice,TNum>().mv_Parameter;

        }

        template<class TLattice, int TNum=1>
        static inline T* getAddress(const int idx) { //Returns pointer to parameter at lattice point k and
                                             //direction 0
            //static Parameter<TObj,TLattice,T,TNumPrefactor*TNum>& instance = getInstance<TLattice,TNum>();
            return &getInstance<TLattice,TNum>().mv_Parameter[idx];

        }

        template<class TLattice, int TNum=1>
        static inline T* getAddress(int idx1, int idx2) { //Returns pointer to parameter at lattice point k and
                                             //direction 0
            //static Parameter<TObj,TLattice,T,TNumPrefactor*TNum>& instance = getInstance<TLattice,TNum>();
            return &getInstance<TLattice,TNum>().mv_Parameter[idx1*instances*TNum + idxlambda<TNum>(idx2)];

        }

        template<class TLattice, int TNum=1>
        static inline T* getAddress(int idx1, int idx2, int idx3) { //Returns pointer to parameter at lattice point k and
                                             //direction 0
            static Parameter<TObj,TLattice,T,TNumPrefactor*TNum>& instance = getInstance<TLattice,TNum>();
            return &instance.mv_Parameter[idx1*instances*TNum + (instances>1)*idx2*TNum + (TNum>1)*idx3];

        }

        template<class TLattice, int TNum=1>
        static inline T& get(const int idx) { //Returns const parameter at index idx

            //static Parameter<TObj,TLattice,T,TNumPrefactor*TNum>& instance = getInstance<TLattice,TNum>();
            return getInstance<TLattice,TNum>().mv_Parameter[idx];
            
        }

        template<class TLattice, int TNum=1>
        static inline T& get(int idx1, int idx2) { //Returns const parameter at index idx
 //MAKE MORE EFFICIENT!!!
 
            //static Parameter<TObj,TLattice,T,TNumPrefactor*TNum>& instance = getInstance<TLattice,TNum>();
            return getInstance<TLattice,TNum>().mv_Parameter[idx1*instances*TNum + idxlambda<TNum>(idx2)];
            
        }

        template<class TLattice, int TNum=1>
        static inline T& get(int idx1, int idx2, int idx3) { //Returns const parameter at index idx
 //MAKE MORE EFFICIENT!!!
            //static Parameter<TObj,TLattice,T,TNumPrefactor*TNum>& instance = getInstance<TLattice,TNum>();
            return getInstance<TLattice,TNum>().mv_Parameter[idx1*instances*TNum + (instances>1)*idx2*TNum + (TNum>1)*idx3];
            
        }

        template<class TLattice, int TNum=1>
        static inline void initialise(const T val,const int idx1, const int idx2, const int idx3){
            #pragma omp critical
            {
            if (!getInstance<TLattice,TNum>().mmInitialised.count(idx1*instances*TNum + (instances>1)*idx2*TNum + (TNum>1)*idx3)) {
                getInstance<TLattice,TNum>().mv_Parameter[idx1*instances*TNum + (instances>1)*idx2*TNum + (TNum>1)*idx3]=val;
                getInstance<TLattice,TNum>().mmInitialised.insert({idx1*instances*TNum + (instances>1)*idx2*TNum + (TNum>1)*idx3, true});
            }
            else {
                
                getInstance<TLattice,TNum>().mmInitialised.erase(idx1*instances*TNum + (instances>1)*idx2*TNum + (TNum>1)*idx3); 
                
            }   
            }        
        }

        template<class TLattice, int TNum=1>
        static inline void initialise(const T val,const int idx1, const int idx2=0){
            #pragma omp critical
            {
            if (!getInstance<TLattice,TNum>().mmInitialised.count(idx1*instances*TNum + idxlambda<TNum>(idx2))) {
                getInstance<TLattice,TNum>().mv_Parameter[idx1*instances*TNum + idxlambda<TNum>(idx2)]=val;
                getInstance<TLattice,TNum>().mmInitialised.insert({idx1*instances*TNum + idxlambda<TNum>(idx2), true});
            }
            else {
                
                getInstance<TLattice,TNum>().mmInitialised.erase(idx1*instances*TNum + idxlambda<TNum>(idx2)); 
                
            }   
            }        
        }

        template<class TLattice, int TNum=1, int idx1=0, int idx2=0>
        static void set(T val) {

            for(int k = TLattice::HaloSize; k < TLattice::N-TLattice::HaloSize; k++) {

                if constexpr(idx2==0) initialise<TLattice,TNum>(val,k,idx1);
                else initialise<TLattice,TNum>(val,k,idx1,idx2);

            }
        }


        template<class TLattice, int TNum=1, int idx1=0, int idx2=0>
        static void set(T (*condition)(const int)) {

            for(int k = TLattice::HaloSize; k < TLattice::N-TLattice::HaloSize; k++) {
                
                if constexpr(idx2==0) initialise<TLattice,TNum>(condition(k),k,idx1);
                else initialise<TLattice,TNum>(condition(k),k,idx1,idx2);

            }

        }



        template<class TLattice, int TNum=1, int idx1=0, int idx2=0>
        static void set(bool (*condition)(const int), T val) {

            for(int k = TLattice::HaloSize; k < TLattice::N-TLattice::HaloSize; k++) {

                if constexpr(idx2==0) {
                    if (condition(k)) initialise<TLattice,TNum>(val,k,idx1);
                }
                else if (condition(k)) initialise<TLattice,TNum>(val,k,idx1,idx2);

            }

        }

        template<class TLattice, int TNum=1, int idx1=0, int idx2=0>
        static void set(bool (*condition)(const int), T val, T false_val) {

            for(int k = TLattice::HaloSize; k < TLattice::N-TLattice::HaloSize; k++) {

                if (condition(k)) {
                    if constexpr(idx2==0) initialise<TLattice,TNum>(val,k,idx1,idx2);
                    else initialise<TLattice,TNum>(false_val,k,idx1,idx2);
                }
                else {
                    if constexpr(idx2==0) initialise<TLattice,TNum>(false_val,k,idx1);
                    else initialise<TLattice,TNum>(false_val,k,idx1,idx2);
                }

            }

        }


        ParameterSingleton(ParameterSingleton<TObj,T,TNumPrefactor> const&)=delete;

        void operator=(ParameterSingleton<TObj,T,TNumPrefactor> const&)=delete;

    private:

        ParameterSingleton(){};

};

//template<class TObj, class TLattice, typename T, int TNum>
//std::vector<T> Parameter<TObj, TLattice, T, TNum>::mv_Parameter; //Must allocate memory for static vector outside of class

//template<class TObj, class TLattice, typename T, int TNum>
//std::map<int,bool> Parameter<TObj, TLattice, T, TNum>::mmInitialised;

template<class TLattice>
class SaveHandler;

template<class TObj,class TLattice,  typename T, int TNum>
inline void Parameter<TObj, TLattice, T, TNum>::save(int t) { //Function to save parameter stored in this class
    SaveHandler<TLattice>::template saveParameter<TObj,TNum>(t);
}

template <int TNDIM>
struct Boundary{
    int Id;
    bool IsCorner;
    std::array<int8_t,TNDIM> NormalDirection;
};

template<int TNDIM,int TInstances = 1>
struct BoundaryLabels : public ParameterSingleton<BoundaryLabels<TNDIM,TInstances>,Boundary<TNDIM>,TInstances> {
    static constexpr const char *mName = "BoundaryLabels";
}; //Labelling of geometry


template<int TInstances=1>
struct Velocity : public ParameterSingleton<Velocity<TInstances>, double, TInstances> {
    
    static constexpr const char *mName = "Velocity";

}; //Velocity, with directions D corresponding to the number of cartesian directions in the stencil

template<int TInstances=1>
struct VelocityOld : public ParameterSingleton<VelocityOld<TInstances>, double, TInstances> {
    
    static constexpr const char *mName = "VelocityOld";

};

template<int TInstances=1>
struct Density : public ParameterSingleton<Density<TInstances>, double, TInstances> {

    static constexpr const char *mName = "Density";

}; //Density

template<int TInstances=1>
struct DensityOld : public ParameterSingleton<DensityOld<TInstances>, double, TInstances> {

    static constexpr const char *mName = "DensityOld";

}; //Density

template<int TInstances=1>
struct Pressure : public ParameterSingleton<Pressure<TInstances>, double, TInstances>{

    static constexpr const char *mName = "Pressure";

}; //Presure

template<int TInstances=1>
struct PressureOld : public ParameterSingleton<PressureOld<TInstances>, double, TInstances>{

    static constexpr const char *mName = "PressureOld";

}; //Presure

template<int TInstances=1>
struct OrderParameter : public ParameterSingleton<OrderParameter<TInstances>, double, TInstances> {

    static constexpr const char *mName = "OrderParameter";

}; //Order parameter representing relative concentration of the phases

template<int TInstances=1>
struct ChemicalPotential : public ParameterSingleton<ChemicalPotential<TInstances>, double, TInstances> {

    static constexpr const char *mName = "ChemicalPotential";

}; //Chemical potential for the multicomponent model

template<int TInstances=1>
struct ChemicalPotentialOld : public ParameterSingleton<ChemicalPotentialOld<TInstances>, double, TInstances> {

    static constexpr const char *mName = "ChemicalPotentialOld";

}; //Chemical potential for the multicomponent model

template<int TInstances=1>
struct OrderParameterOld : public ParameterSingleton<OrderParameterOld<TInstances>, double, TInstances> {

    static constexpr const char *mName = "OrderParameterOld";

};

template<int TInstances=1>
struct Humidity : public ParameterSingleton<Humidity<TInstances>> {

    static constexpr const char *mName = "Humidity";

};

template<int TInstances=1>
struct Temperature : public ParameterSingleton<Temperature<TInstances>> {

    static constexpr const char *mName = "Humidity";

};

template<int TInstances=1>
struct HumidityOld : public ParameterSingleton<HumidityOld<TInstances>> {

    static constexpr const char *mName = "HumidityOld";

};

template<int TInstances=1>
struct MassSink : public ParameterSingleton<MassSink<TInstances>> {

    static constexpr const char *mName = "MassSink";

};

template<class TObj, int TInstances = 1>
struct Laplacian : public ParameterSingleton<Laplacian<TObj,TInstances>,double,TInstances> {
    static constexpr char mName[9+sizeof(TObj::mName)] = "Laplacian"+TObj::mName;
}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct LaplacianChemicalPotential : public Laplacian<ChemicalPotential<TInstances>,TInstances> {

    static constexpr const char *mName = "LaplacianChemicalPotential";

}; //Laplacian of the order parameter

template<int TInstances = 1>
struct LaplacianDensity : public Laplacian<Density<TInstances>,TInstances> {

    static constexpr const char *mName = "LaplacianDensity";

}; //Laplacian of the order parameter

template<int TInstances = 1>
struct LaplacianOrderParameter : public Laplacian<OrderParameter<TInstances>,TInstances> {

    static constexpr const char *mName = "LaplacianOrderParameter";

}; //Laplacian of the order parameter

template<class TObj, int TInstances=1>
struct Gradient : public ParameterSingleton<Gradient<TObj,TInstances>, double, TInstances> {
    static constexpr char mName[8+sizeof(TObj::mName)] = "Gradient"+TObj::mName;
}; //Directional first order gradients of the order parameter

template<class TObj, int TInstances = 1>
struct GradientMixed : public ParameterSingleton<GradientMixed<TObj,TInstances>, double, TInstances> {
    static constexpr char mName[13+sizeof(TObj::mName)] = "GradientMixed"+TObj::mName;
}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct GradientOrderParameter : public Gradient<OrderParameter<TInstances>,TInstances> {

    static constexpr const char *mName="GradientOrderParameter";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct GradientDensity : public Gradient<Density<TInstances>,TInstances> {

    static constexpr const char *mName="GradientDensity";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct GradientChemicalPotential : public Gradient<ChemicalPotential<TInstances>,TInstances> {

    static constexpr const char *mName="GradientChemicalPotential";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct GradientPressure : public Gradient<Pressure<TInstances>,TInstances> {

    static constexpr const char *mName="GradientPressure";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct GradientHumidity : public Gradient<Humidity<TInstances>,TInstances> {

    static constexpr const char *mName="GradientHumidity";

};

template<int TInstances = 1>
struct MixedGradientOrderParameter : public GradientMixed<OrderParameter<TInstances>,TInstances> {

    static constexpr const char *mName="MixedGradientOrderParameter";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct MixedGradientDensity : public GradientMixed<Density<TInstances>,TInstances> {

    static constexpr const char *mName="MixedGradientDensity";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct MixedGradientPressure : public GradientMixed<Pressure<TInstances>,TInstances> {

    static constexpr const char *mName="MixedGradientPressure";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct Tau : public ParameterSingleton<Tau<TInstances>,double,TInstances> {
    static constexpr const char *mName = "Tau";
}; //Labelling of geometry

template<int TInstances = 1>
struct InverseTau : public ParameterSingleton<InverseTau<TInstances>,double,TInstances> {
    static constexpr const char *mName = "InverseTau";
}; //Labelling of geometry

template<int TInstances = 1>
struct LogFugacity : public ParameterSingleton<LogFugacity<TInstances>,double,TInstances> {
    static constexpr const char *mName = "Fugacity";
}; //Labelling of geometry

template<int TInstances = 1>
struct GradientLogFugacity : public ParameterSingleton<GradientLogFugacity<TInstances>,double,TInstances> {
    static constexpr const char *mName = "Fugacity";
}; //Labelling of geometry
