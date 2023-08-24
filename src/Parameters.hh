#pragma once
#include <vector>
#include<string>
#include<map>
#include<any>
#include "Global.hh"
#include <fstream>
#include <iostream>

//Parameters.hh: This file details how macroscopic quantities are stored and interacted with.
//The Distribution class contains some vectors and functions to return values from these. This will be inherited
//from in the Data classes (see Data.hh). The Parameter class contains a static vector and functions to return
//the parameters stored. The reason for this static vector is to ensure that, if I use "Density" in multiple
//classes, the value will be consistent between classes. Note that, for every template configuration, a new
//static variable is created, so I pass each class to itself as a template to ensure the parameters are unique.

template<class TStencil> //Distribution must know information about the stencil as this determines the size
                         //of the vectors and the data layout
struct Distribution_Base { //Distribution base class

    Distribution_Base(std::vector<int>& neighbors) : mv_DistNeighbors(neighbors) {

        for(int idx = 0; idx <TStencil::Q; idx++) { //Calculate the k offset for the neighbors in each direction

            ma_Opposites[idx] =TStencil::Opposites[idx];

        }

    }

    inline void saveEquilibriums(const double* equilibrium, int k) {}

    /**
    * \brief Returns the opposite index at the chosen index (Rotation by 180 degrees).
    */

    int ma_Opposites[TStencil::Q]; //!<Array containing the opposite indices at each index (rotated by 180 degrees).

    inline int getOpposite(int idx) {

        return ma_Opposites[idx];

    }

    std::vector<double> mv_Distribution; //Vector that will store the distribution information
    std::vector<double> mv_OldDistribution; //Possibly unused vector storing the old distributions (at t_old=t_new-1)
                                       //This is required in the collision step
    std::vector<double> mv_CommDistribution; //Vector that will store the distribution information for communication
    std::vector<double> mv_EquilibriumDistribution;
    std::vector<int>& mv_DistNeighbors; //Reference to vector containing neighbor information
    inline std::vector<double>& getDistribution() { //Get a vector containing the distributions

        return mv_Distribution;

    }

    inline std::vector<double>& getCommDistribution() { //Get a vector containing the distributions for communication (along direction `neighbor`?)
        return mv_CommDistribution;
    }

    inline const double* getDistributionPointer(const int k) const { //Get a constant pointer to the the distribution at
                                                             //lattice point k and pointing in direction 0
        return &mv_Distribution[k * TStencil::Q];
    }
    inline double* getDistributionPointer(const int k) { //Get a pointer to the the distribution at
                                                 //lattice point k and pointing in direction 0
        return &mv_Distribution[k * TStencil::Q];

    }
    inline const double& getDistribution(const int idx) const { //Get const distribution value at a given index

        return mv_Distribution[idx];

    }
    inline double& getDistribution(const int idx) { //Get distribution value at a given index

        return mv_Distribution[idx];

    }
    inline std::vector<double>& getDistributionOld() { //SAME AS ABOVE BUT FOR OLD DISTRIBUTION

        return mv_OldDistribution;

    }
    inline const double* getDistributionOldPointer(const int k) const {

        return &mv_OldDistribution[k *TStencil::Q];

    }
    inline double* getDistributionOldPointer(const int k) {

        return &mv_OldDistribution[k *TStencil::Q];

    }
    inline const double& getDistributionOld(const int k) const {

        return mv_OldDistribution[k];

    }
    inline double& getDistributionOld(const int k) {

        return mv_OldDistribution[k];

    }
    inline std::vector<double>& getEquilibrium() { //Get a vector containing the distributions

        return mv_EquilibriumDistribution;

    }
    inline const double* getEquilibriumPointer(const int k) const { //Get a constant pointer to the the distribution at
                                                             //lattice point k and pointing in direction 0
        return &mv_EquilibriumDistribution[k *TStencil::Q];
    }
    inline double* getEquilibriumPointer(const int k) { //Get a pointer to the the distribution at
                                                 //lattice point k and pointing in direction 0
        return &mv_EquilibriumDistribution[k *TStencil::Q];

    }
    inline const double& getEquilibrium(const int idx) const { //Get const distribution value at a given index

        return mv_EquilibriumDistribution[idx];

    }
    inline double& getEquilibrium(const int idx) { //Get distribution value at a given index

        return mv_EquilibriumDistribution[idx];

    }

    int mQ =TStencil::Q;

    using Stencil =TStencil;
};

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

        using ParamType = T;

        inline void Save(std::string filename, int t, std::string datadir);

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

            return getInstance<TLattice,TNum>().mv_Parameter;

        }

        template<class TLattice, int TNum=1>
        static inline T* getAddress(const int idx) { //Returns pointer to parameter at lattice point k and
                                             //direction 0
            return &getInstance<TLattice,TNum>().mv_Parameter[idx];

        }

        template<class TLattice, int TNum=1>
        static inline T* getAddress(int idx1, int idx2) { //Returns pointer to parameter at lattice point k and
                                             //direction 0
            return &getInstance<TLattice,TNum>().mv_Parameter[idx1*instances*TNum + idxlambda<TNum>(idx2)];

        }

        template<class TLattice, int TNum=1>
        static inline T* getAddress(int idx1, int idx2, int idx3) { //Returns pointer to parameter at lattice point k and
                                             //direction 0
            return &getInstance<TLattice,TNum>().mv_Parameter[idx1*instances*TNum + (instances>1)*idx2*TNum + (TNum>1)*idx3];

        }

        template<class TLattice, int TNum=1>
        static inline T& get(const int idx) { //Returns const parameter at index idx

            return getInstance<TLattice,TNum>().mv_Parameter[idx];
            
        }

        template<class TLattice, int TNum=1>
        static inline T& get(int idx1, int idx2) { //Returns const parameter at index idx
 //MAKE MORE EFFICIENT!!!
 
            return getInstance<TLattice,TNum>().mv_Parameter[idx1*instances*TNum + idxlambda<TNum>(idx2)];
            
        }

        template<class TLattice, int TNum=1>
        static inline T& get(int idx1, int idx2, int idx3) { //Returns const parameter at index idx
 //MAKE MORE EFFICIENT!!!
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

template<class TObj,class TLattice,  typename T, int TNum>
inline void Parameter<TObj, TLattice, T, TNum>::Save(std::string filename, int t, std::string datadir) { //Function to save parameter stored in this class

    char fdump[512];
    sprintf(fdump, "%s/%s_t%i.mat", datadir.c_str(), filename.c_str(), t); //Buffer containing file name and location.

#ifdef MPIPARALLEL //When MPI is used we need a different approach for saving as all nodes are trying to write to the file

    MPI_File fh;

    MPI_File_open(MPI_COMM_SELF, fdump, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh); //Open the file using mpi in write only mode
    
    MPI_File_seek(fh, sizeof(typename TObj::ParamType) * mpi.rank * mNum * (TLattice::LX * TLattice::LY * TLattice::LZ) / mpi.size, MPI_SEEK_SET); //Skip to a certain location in the file, currently
    
    MPI_File_write(fh,&mv_Parameter[TLattice::HaloSize * mNum], mNum * (TLattice::N - 2 * TLattice::HaloSize), MPI_DOUBLE, MPI_STATUSES_IGNORE);

    MPI_File_close(&fh);
         
#else

    std::ofstream fs(fdump, std::ios::out | std::ios::binary);
    
    fs.seekp(sizeof(typename TObj::ParamType) * mpi.rank * mNum * (TLattice::LX * TLattice::LY * TLattice::LZ) / mpi.size);

    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) { 

        for(int idx = 0; idx < mNum; idx++) {
            //std::cout<<*(char *)(&mv_Parameter[k * mNum + idx])<<std::endl;
            fs.write((char *)(&mv_Parameter[k * mNum + idx]), sizeof(typename TObj::ParamType));
        }
        
    };

    fs.close();

#endif

}

struct Boundary{
    int Id;
    bool IsCorner;
    std::vector<int8_t> NormalDirection;
};

template<int TInstances = 1>
struct BoundaryLabels : public ParameterSingleton<BoundaryLabels<TInstances>,Boundary,TInstances> {
    static constexpr char mName[] = "BoundaryLabels";
}; //Labelling of geometry

template<class TLattice>
class ParameterSave {

    public:

        ParameterSave(std::string datadir) : mDataDir(datadir){

            int status = system(((std::string)"mkdir -p " + mDataDir).c_str());
            if (status) std::cout << "Error creating output directory" << std::endl;

        #ifdef MPIPARALLEL

            MPI_Type_create_resized(MPI_INT, 0L, sizeof(Boundary), &mMPIBoundary)
            MPI_Type_commit(&mMPIBoundary);

        #endif

        }
        ParameterSave(ParameterSave& other) : mDataDir(other.datadir){

            int status = system(((std::string)"mkdir -p " + mDataDir).c_str());
            if (status) std::cout << "Error creating output directory" << std::endl;

        #ifdef MPIPARALLEL

            MPI_Type_create_resized(MPI_INT, 0L, sizeof(Boundary), &mMPIBoundary)
            MPI_Type_commit(&mMPIBoundary);

        #endif

        }

        inline void SaveHeader(const int& timestep, const int& saveinterval);

        template<class TParameter, int TNumDir=1>
        void SaveParameter(int timestep){
            TParameter::template getInstance<TLattice,TNumDir>().Save(TParameter::mName,timestep,mDataDir);
        }

        template<class... TParameter>
        void SaveParameter(int timestep, TParameter&... params){
            (params.Save(TParameter::mName,timestep,mDataDir),...);
        }

        void SaveBoundaries(int timestep){
            char fdump[512];
            sprintf(fdump, "%s/%s_t%i.mat", mDataDir.c_str(), BoundaryLabels<>::mName, timestep); //Buffer containing file name and location.

        #ifdef MPIPARALLEL //When MPI is used we need a different approach for saving as all nodes are trying to write to the file

            MPI_File fh;

            MPI_File_open(MPI_COMM_SELF, fdump, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh); //Open the file using mpi in write only mode
            
            MPI_File_seek(fh, sizeof(int) * mpi.rank * (TLattice::LX * TLattice::LY * TLattice::LZ) / mpi.size, MPI_SEEK_SET); //Skip to a certain location in the file, currently
            
            MPI_File_write(fh,&BoundaryLabels<>::get<TLattice>()[TLattice::HaloSize], (TLattice::N - 2 * TLattice::HaloSize), mMPIBoundary, MPI_STATUSES_IGNORE);

            MPI_File_close(&fh);
                
        #else

            std::ofstream fs(fdump, std::ios::out | std::ios::binary);
            
            fs.seekp(sizeof(int) * mpi.rank * (TLattice::LX * TLattice::LY * TLattice::LZ) / mpi.size);

            for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) { 

                
                fs.write((char *)(&BoundaryLabels<>::get<TLattice>()[k].Id), sizeof(int));
                
                
            };

            fs.close();

        #endif

        }

    private:

        std::string mDataDir;

    #ifdef MPIPARALLEL

        MPI_Datatype mMPIBoundary;

    #endif
        
};

template<class TLattice>
inline void ParameterSave<TLattice>::SaveHeader(const int& timestep, const int& saveinterval) { //Function to save parameter stored in this class
    
    if(mpi.rank==0){
        std::cout<<"SAVING HEADER"<<std::endl;
        char fdump[256];
        sprintf(fdump, "%s/Header.mat", mDataDir.c_str()); //Buffer containing file name and location.

        std::ofstream fs(fdump, std::ios::out | std::ios::binary);

        fs.write((char *)(&TLattice::LX), sizeof(int));
        fs.write((char *)(&TLattice::LY), sizeof(int));
        fs.write((char *)(&TLattice::LZ), sizeof(int));
        fs.write((char *)(&TLattice::NDIM), sizeof(int));
        fs.write((char *)(&timestep), sizeof(int));
        fs.write((char *)(&saveinterval), sizeof(int));

        fs.close();

    }

}



template<int TInstances=1>
struct Velocity : public ParameterSingleton<Velocity<TInstances>, double, TInstances> {
    
    static constexpr char mName[] = "Velocity";

}; //Velocity, with directions D corresponding to the number of cartesian directions in the stencil

template<int TInstances=1>
struct VelocityOld : public ParameterSingleton<VelocityOld<TInstances>, double, TInstances> {
    
    static constexpr char mName[] = "VelocityOld";

};

template<int TInstances=1>
struct Density : public ParameterSingleton<Density<TInstances>, double, TInstances> {

    static constexpr char mName[] = "Density";

}; //Density

template<int TInstances=1>
struct Pressure : public ParameterSingleton<Pressure<TInstances>, double, TInstances>{

    static constexpr char mName[] = "Pressure";

}; //Presure

template<int TInstances=1>
struct OrderParameter : public ParameterSingleton<OrderParameter<TInstances>, double, TInstances> {

    static constexpr char mName[] = "OrderParameter";

}; //Order parameter representing relative concentration of the phases

template<int TInstances=1>
struct ChemicalPotential : public ParameterSingleton<ChemicalPotential<TInstances>, double, TInstances> {

    static constexpr char mName[] = "ChemicalPotential";

}; //Chemical potential for the multicomponent model

template<int TInstances=1>
struct OrderParameterOld : public ParameterSingleton<OrderParameterOld<TInstances>, double, TInstances> {

    static constexpr char mName[] = "OrderParameterOld";

};

template<int TInstances=1>
struct Humidity : public ParameterSingleton<Humidity<TInstances>> {

    static constexpr char mName[] = "Humidity";

};

template<int TInstances=1>
struct HumidityOld : public ParameterSingleton<HumidityOld<TInstances>> {

    static constexpr char mName[] = "HumidityOld";

};

template<int TInstances=1>
struct MassSink : public ParameterSingleton<MassSink<TInstances>> {

    static constexpr char mName[] = "MassSink";

};

template<class TObj, int TInstances = 1>
struct Laplacian : public ParameterSingleton<Laplacian<TObj,TInstances>,double,TInstances> {
    static constexpr char mName[9+sizeof(TObj::mName)] = "Laplacian"+TObj::mName;
}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct LaplacianChemicalPotential : public Laplacian<ChemicalPotential<TInstances>,TInstances> {

    static constexpr char mName[] = "LaplacianChemicalPotential";

}; //Laplacian of the order parameter

template<int TInstances = 1>
struct LaplacianDensity : public Laplacian<Density<TInstances>,TInstances> {

    static constexpr char mName[] = "LaplacianDensity";

}; //Laplacian of the order parameter

template<int TInstances = 1>
struct LaplacianOrderParameter : public Laplacian<OrderParameter<TInstances>,TInstances> {

    static constexpr char mName[] = "LaplacianOrderParameter";

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

    static constexpr char mName[]="GradientOrderParameter";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct GradientDensity : public Gradient<Density<TInstances>,TInstances> {

    static constexpr char mName[]="GradientDensity";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct GradientChemicalPotential : public Gradient<ChemicalPotential<TInstances>,TInstances> {

    static constexpr char mName[]="GradientChemicalPotential";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct GradientPressure : public Gradient<Pressure<TInstances>,TInstances> {

    static constexpr char mName[]="GradientPressure";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct GradientHumidity : public Gradient<Humidity<TInstances>,TInstances> {

    static constexpr char mName[]="GradientHumidity";

};

template<int TInstances = 1>
struct MixedGradientOrderParameter : public GradientMixed<OrderParameter<TInstances>,TInstances> {

    static constexpr char mName[]="MixedGradientOrderParameter";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct MixedGradientDensity : public GradientMixed<Density<TInstances>,TInstances> {

    static constexpr char mName[]="MixedGradientDensity";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct MixedGradientPressure : public GradientMixed<Pressure<TInstances>,TInstances> {

    static constexpr char mName[]="MixedGradientPressure";

}; //Directional first order gradients of the order parameter

template<int TInstances = 1>
struct Tau : public ParameterSingleton<Tau<TInstances>,double,TInstances> {
    static constexpr char mName[] = "Tau";
}; //Labelling of geometry

template<int TInstances = 1>
struct InverseTau : public ParameterSingleton<InverseTau<TInstances>,double,TInstances> {
    static constexpr char mName[] = "InverseTau";
}; //Labelling of geometry

template<int TInstances = 1>
struct LogFugacity : public ParameterSingleton<LogFugacity<TInstances>,double,TInstances> {
    static constexpr char mName[] = "Fugacity";
}; //Labelling of geometry

template<int TInstances = 1>
struct GradientLogFugacity : public ParameterSingleton<GradientLogFugacity<TInstances>,double,TInstances> {
    static constexpr char mName[] = "Fugacity";
}; //Labelling of geometry
