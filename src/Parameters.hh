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

template<class T_stencil> //Distribution must know information about the stencil as this determines the size
                         //of the vectors and the data layout
struct Distribution_Base { //Distribution base class

    Distribution_Base(std::vector<int>& neighbors) : mv_DistNeighbors(neighbors) {

        for(int idx = 0; idx < T_stencil::Q; idx++) { //Calculate the k offset for the neighbors in each direction

            ma_Opposites[idx] = T_stencil::Opposites[idx];

        }

    }

    /**
    * \brief Returns the opposite index at the chosen index (Rotation by 180 degrees).
    */

    int ma_Opposites[T_stencil::Q]; //!<Array containing the opposite indices at each index (rotated by 180 degrees).

    inline int getOpposite(int idx) {

        return ma_Opposites[idx];

    }

    std::vector<double> mv_Distribution; //Vector that will store the distribution information
    std::vector<double> mv_OldDistribution; //Possibly unused vector storing the old distributions (at t_old=t_new-1)
                                       //This is required in the collision step
    std::vector<int>& mv_DistNeighbors; //Reference to vector containing neighbor information
    inline std::vector<double>& getDistribution() { //Get a vector containing the distributions

        return mv_Distribution;

    }
    inline const double* getDistributionPointer(const int k) const { //Get a constant pointer to the the distribution at
                                                             //lattice point k and pointing in direction 0
        return &mv_Distribution[k * T_stencil::Q];
    }
    inline double* getDistributionPointer(const int k) { //Get a pointer to the the distribution at
                                                 //lattice point k and pointing in direction 0
        return &mv_Distribution[k * T_stencil::Q];

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

        return &mv_OldDistribution[k * T_stencil::Q];

    }
    inline double* getDistributionOldPointer(const int k) {

        return &mv_OldDistribution[k * T_stencil::Q];

    }
    inline const double& getDistributionOld(const int k) const {

        return mv_OldDistribution[k];

    }
    inline double& getDistributionOld(const int k) {

        return mv_OldDistribution[k];

    }

    int m_Q = T_stencil::Q;

    using Stencil = T_stencil;
};

template<class T_obj, class T_lattice, typename T, int num=1> //obj template will guarantee a unique instance of the class with its own
                                       //static vector. I pass the class to itself to guarantee this
                                       //
                                       //T determines the type stored in the vector, num
                                       //determines the number of directions for each lattice point (eg you might
                                       //want two directions corresponding to velocity in x and velocity in y at
                                       //each point.
class Parameter {

    public:

        std::map<int,bool> mm_Initialised; //Change this to std::set

        static constexpr int m_Num=num;
        static constexpr int m_Instances=T_obj::instances;
        static constexpr int m_Directions=m_Num/T_obj::instances;

        using ParamType = T;

        inline void Save(std::string filename, int t, std::string datadir);

        std::vector<T> mv_Parameter; //Static vector (Does not change between objects of the class)

        template<class, class, int>
        friend class ParameterSingleton;

    private:

        Parameter() {

            mv_Parameter.resize(num * T_lattice::N); //Resize to the desired size

        }
        
        Parameter(const Parameter& other) {}

        
        
};

template<class T_obj, typename T=double, int numprefactor=1>
class ParameterSingleton {

    public:
        static constexpr int instances=numprefactor;
        static constexpr bool multipleinstances = (instances>1);
        
        template<class T_lattice, int num=1>
        static inline Parameter<T_obj,T_lattice,T,numprefactor*num>& getInstance() { static Parameter<T_obj,T_lattice,T,numprefactor*num> instance;
                                                                          return instance;};

        template<class T_lattice, int num=1>
        static inline std::vector<T>& get() { //Returns vector containing the parameter

            return getInstance<T_lattice,num>().mv_Parameter;

        }

        template<class T_lattice, int num=1>
        static inline T* getAddress(const int idx) { //Returns pointer to parameter at lattice point k and
                                             //direction 0
            return &getInstance<T_lattice,num>().mv_Parameter[idx];

        }

        template<class T_lattice, int num=1>
        static inline T* getAddress(int idx1, int idx2) { //Returns pointer to parameter at lattice point k and
                                             //direction 0
            constexpr bool multipledirections = (num>1);

            return &getInstance<T_lattice,num>().mv_Parameter[idx1*instances*num+multipleinstances*(idx2)+(!multipleinstances)*multipledirections*idx2];

        }

        template<class T_lattice, int num=1>
        static inline T* getAddress(int idx1, int idx2, int idx3) { //Returns pointer to parameter at lattice point k and
                                             //direction 0
            constexpr bool multipledirections = (num>1);

            return &getInstance<T_lattice,num>().mv_Parameter[idx1*instances*num+multipleinstances*(idx2*num+multipledirections*idx3)+(!multipleinstances)*multipledirections*idx3];

        }

        template<class T_lattice, int num=1>
        static inline T& get(const int idx) { //Returns const parameter at index idx

            return getInstance<T_lattice,num>().mv_Parameter[idx];
            
        }

        template<class T_lattice, int num=1>
        static inline T& get(int idx1, int idx2) { //Returns const parameter at index idx

            constexpr bool multipledirections = (num>1);
 //MAKE MORE EFFICIENT!!!
            return getInstance<T_lattice,num>().mv_Parameter[idx1*instances*num+multipleinstances*(idx2)+(!multipleinstances)*multipledirections*idx2];
            
        }

        template<class T_lattice, int num=1>
        static inline T& get(int idx1, int idx2, int idx3) { //Returns const parameter at index idx

            constexpr bool multipledirections = (num>1);
 //MAKE MORE EFFICIENT!!!
            return getInstance<T_lattice,num>().mv_Parameter[idx1*instances*num+multipleinstances*(idx2*num+multipledirections*idx3)+(!multipleinstances)*multipledirections*idx3];
            
        }

        template<class T_lattice, int num=1>
        static inline void initialise(const T val,const int idx1, const int idx2=0, const int idx3=0){
            constexpr bool multipledirections = (num>1);
            #pragma omp critical
            {
            if (!getInstance<T_lattice,num>().mm_Initialised.count(idx1*instances*num+multipleinstances*(idx2*num+multipledirections*idx3)+(!multipleinstances)*multipledirections*idx3)) {
                getInstance<T_lattice,num>().mv_Parameter[idx1*instances*num+multipleinstances*(idx2*num+multipledirections*idx3)+(!multipleinstances)*multipledirections*idx3]=val;
                getInstance<T_lattice,num>().mm_Initialised.insert({idx1*instances*num+multipleinstances*(idx2*num+multipledirections*idx3)+(!multipleinstances)*multipledirections*idx3,true});
            }
            else {
                
                getInstance<T_lattice,num>().mm_Initialised.erase(idx1*instances*num+multipleinstances*(idx2*num+multipledirections*idx3)+(!multipleinstances)*multipledirections*idx3); 
                
            }   
            }        
        }

        template<class T_lattice, int num=1, int idx1=0, int idx2=0>
        static void set(T (*condition)(const int)) {

            for(int k = T_lattice::HaloSize; k < T_lattice::N-T_lattice::HaloSize; k++) {

                initialise<T_lattice,num>(condition(k),k,idx1,idx2);

            }

        }

        template<class T_lattice, int num=1, int idx1=0, int idx2=0>
        static void set(bool (*condition)(const int), T val) {

            for(int k = T_lattice::HaloSize; k < T_lattice::N-T_lattice::HaloSize; k++) {

                if (condition(k)) initialise<T_lattice,num>(val,k,idx1,idx2);

            }

        }

        template<class T_lattice, int num=1, int idx1=0, int idx2=0>
        static void set(bool (*condition)(const int), T val, T false_val) {

            for(int k = T_lattice::HaloSize; k < T_lattice::N-T_lattice::HaloSize; k++) {

                if (condition(k)) initialise<T_lattice,num>(val,k,idx1,idx2);
                else initialise<T_lattice,num>(false_val,k,idx1,idx2);

            }

        }


        ParameterSingleton(ParameterSingleton<T_obj,T,numprefactor> const&)=delete;

        void operator=(ParameterSingleton<T_obj,T,numprefactor> const&)=delete;

    private:

        ParameterSingleton(){};


};

//template<class T_obj, class T_lattice, typename T, int num>
//std::vector<T> Parameter<T_obj, T_lattice, T, num>::mv_Parameter; //Must allocate memory for static vector outside of class

//template<class T_obj, class T_lattice, typename T, int num>
//std::map<int,bool> Parameter<T_obj, T_lattice, T, num>::mm_Initialised;

template<class T_obj,class T_lattice,  typename T, int num>
inline void Parameter<T_obj, T_lattice, T, num>::Save(std::string filename, int t, std::string datadir) { //Function to save parameter stored in this class

    char fdump[512];
    sprintf(fdump, "%s/%s_t%i.mat", datadir.c_str(), filename.c_str(), t); //Buffer containing file name and location.

#ifdef MPIPARALLEL //When MPI is used we need a different approach for saving as all nodes are trying to write to the file

    MPI_File fh;

    MPI_File_open(MPI_COMM_SELF, fdump, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh); //Open the file using mpi in write only mode
    
    MPI_File_seek(fh, sizeof(double) * mpi.rank * m_Num * (T_lattice::LX * T_lattice::LY * T_lattice::LZ) / mpi.size, MPI_SEEK_SET); //Skip to a certain location in the file, currently
    
    MPI_File_write(fh,&mv_Parameter[T_lattice::HaloSize * m_Num], m_Num * (T_lattice::N - 2 * T_lattice::HaloSize), MPI_DOUBLE, MPI_STATUSES_IGNORE);

    MPI_File_close(&fh);
         
#else

    std::ofstream fs(fdump, std::ios::out | std::ios::binary);
    
    fs.seekp(sizeof(double) * mpi.rank * m_Num * (T_lattice::LX * T_lattice::LY * T_lattice::LZ) / mpi.size);

    for (int k = T_lattice::HaloSize; k < T_lattice::N - T_lattice::HaloSize; k++) { 

        for(int idx = 0; idx < m_Num; idx++) fs.write((char *)(&mv_Parameter[k * m_Num + idx]), sizeof(double));
        
    };

    fs.close();

#endif

}

template<class T_lattice>
class ParameterSave {

    public:

        ParameterSave(std::string datadir) : m_DataDir(datadir){

            int status = system(((std::string)"mkdir -p " + m_DataDir).c_str());
            if (status) std::cout << "Error creating output directory" << std::endl;

        }
        ParameterSave(ParameterSave& other) : m_DataDir(other.datadir){

            int status = system(((std::string)"mkdir -p " + m_DataDir).c_str());
            if (status) std::cout << "Error creating output directory" << std::endl;

        }

        inline void SaveHeader(const int& timestep, const int& saveinterval);

        template<class T_parameter, int numdir=1>
        void SaveParameter(int timestep){
            T_parameter::template getInstance<T_lattice,numdir>().Save(T_parameter::m_Name,timestep,m_DataDir);
        }

        template<class... T_parameter>
        void SaveParameter(int timestep, T_parameter&... params){
            (params.Save(T_parameter::m_Name,timestep,m_DataDir),...);
        }

    private:

        std::string m_DataDir;
        
};

template<class T_lattice>
inline void ParameterSave<T_lattice>::SaveHeader(const int& timestep, const int& saveinterval) { //Function to save parameter stored in this class
    
    if(mpi.rank==0){
        std::cout<<"SAVING HEADER"<<std::endl;
        char fdump[256];
        sprintf(fdump, "%s/Header.mat", m_DataDir.c_str()); //Buffer containing file name and location.

        std::ofstream fs(fdump, std::ios::out | std::ios::binary);

        fs.write((char *)(&T_lattice::LX), sizeof(int));
        fs.write((char *)(&T_lattice::LY), sizeof(int));
        fs.write((char *)(&T_lattice::LZ), sizeof(int));
        fs.write((char *)(&T_lattice::NDIM), sizeof(int));
        fs.write((char *)(&timestep), sizeof(int));
        fs.write((char *)(&saveinterval), sizeof(int));

        fs.close();

    }

}



template<int instances=1>
struct Velocity : public ParameterSingleton<Velocity<instances>, double, instances> {
    
    static constexpr char m_Name[] = "Velocity";

}; //Velocity, with directions D corresponding to the number of cartesian directions in the stencilUw

template<int instances=1>
struct Density : public ParameterSingleton<Density<instances>, double, instances> {

    static constexpr char m_Name[] = "Density";

}; //Density

template<int instances=1>
struct Pressure : public ParameterSingleton<Pressure<instances>, double, instances>{

    static constexpr char m_Name[] = "Pressure";

}; //Presure

template<int instances=1>
struct OrderParameter : public ParameterSingleton<OrderParameter<instances>, double, instances> {

    static constexpr char m_Name[] = "OrderParameter";

}; //Order parameter representing relative concentration of the phases

template<int instances=1>
struct ChemicalPotential : public ParameterSingleton<ChemicalPotential<instances>> {

    static constexpr char m_Name[] = "ChemicalPotential";

}; //Chemical potential for the multicomponent model

template<class T_obj,int instances = 1>
struct Laplacian : public ParameterSingleton<Laplacian<T_obj,instances>,double,instances> {
    static constexpr char m_Name[9+sizeof(T_obj::m_Name)] = "Laplacian"+T_obj::m_Name;
}; //Directional first order gradients of the order parameter

template<int instances = 1>
struct LaplacianChemicalPotential : public Laplacian<ChemicalPotential<instances>,instances> {

    static constexpr char m_Name[] = "LaplacianChemicalPotential";

}; //Laplacian of the order parameter

template<int instances = 1>
struct LaplacianDensity : public Laplacian<Density<instances>,instances> {

    static constexpr char m_Name[] = "LaplacianDensity";

}; //Laplacian of the order parameter

template<int instances = 1>
struct LaplacianOrderParameter : public Laplacian<OrderParameter<instances>,instances> {

    static constexpr char m_Name[] = "LaplacianOrderParameter";

}; //Laplacian of the order parameter

template<class T_obj,int instances=1>
struct Gradient : public ParameterSingleton<Gradient<T_obj,instances>, double, instances> {
    static constexpr char m_Name[8+sizeof(T_obj::m_Name)] = "Gradient"+T_obj::m_Name;
}; //Directional first order gradients of the order parameter

template<class T_obj,int instances = 1>
struct GradientMixed : public ParameterSingleton<GradientMixed<T_obj,instances>, double, instances> {
    static constexpr char m_Name[13+sizeof(T_obj::m_Name)] = "GradientMixed"+T_obj::m_Name;
}; //Directional first order gradients of the order parameter

template<int instances = 1>
struct GradientOrderParameter : public Gradient<OrderParameter<instances>,instances> {

    static constexpr char m_Name[]="GradientOrderParameter";

}; //Directional first order gradients of the order parameter

template<int instances = 1>
struct GradientDensity : public Gradient<Density<instances>,instances> {

    static constexpr char m_Name[]="GradientDensity";

}; //Directional first order gradients of the order parameter

template<int instances = 1>
struct GradientChemicalPotential : public Gradient<ChemicalPotential<instances>,instances> {

    static constexpr char m_Name[]="GradientChemicalPotential";

}; //Directional first order gradients of the order parameter

template<int instances = 1>
struct GradientPressure : public Gradient<Pressure<instances>,instances> {

    static constexpr char m_Name[]="GradientPressure";

}; //Directional first order gradients of the order parameter

template<int instances = 1>
struct MixedGradientOrderParameter : public GradientMixed<OrderParameter<instances>,instances> {

    static constexpr char m_Name[]="MixedGradientOrderParameter";

}; //Directional first order gradients of the order parameter

template<int instances = 1>
struct MixedGradientDensity : public GradientMixed<Density<instances>,instances> {

    static constexpr char m_Name[]="MixedGradientDensity";

}; //Directional first order gradients of the order parameter

template<int instances = 1>
struct MixedGradientPressure : public GradientMixed<Pressure<instances>,instances> {

    static constexpr char m_Name[]="MixedGradientPressure";

}; //Directional first order gradients of the order parameter

template<int instances = 1>
struct SolidLabels : public ParameterSingleton<SolidLabels<instances>,int,instances> {
    static constexpr char m_Name[] = "SolidLabels";
}; //Labelling of geometry

template<int instances = 1>
struct Tau : public ParameterSingleton<Tau<instances>,double,instances> {
    static constexpr char m_Name[] = "Tau";
}; //Labelling of geometry

template<int instances = 1>
struct InverseTau : public ParameterSingleton<InverseTau<instances>,double,instances> {
    static constexpr char m_Name[] = "InverseTau";
}; //Labelling of geometry

template<int instances = 1>
struct LogFugacity : public ParameterSingleton<LogFugacity<instances>,double,instances> {
    static constexpr char m_Name[] = "Fugacity";
}; //Labelling of geometry

template<int instances = 1>
struct GradientLogFugacity : public ParameterSingleton<GradientLogFugacity<instances>,double,instances> {
    static constexpr char m_Name[] = "Fugacity";
}; //Labelling of geometry
