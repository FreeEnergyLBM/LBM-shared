#pragma once
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

template<class stencil> //Distribution must know information about the stencil as this determines the size
                         //of the vectors and the data layout
struct Distribution_Base { //Distribution base class

    Distribution_Base(std::vector<int>& neighbors) : mv_DistNeighbors(neighbors) {

        for(int idx = 0; idx < stencil::Q; idx++) { //Calculate the k offset for the neighbors in each direction

            ma_Opposites[idx] = stencil::Opposites[idx];

        }

    }

    /**
    * \brief Returns the opposite index at the chosen index (Rotation by 180 degrees).
    */

    int ma_Opposites[stencil::Q]; //!<Array containing the opposite indices at each index (rotated by 180 degrees).

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
        return &mv_Distribution[k * stencil::Q];
    }
    inline double* getDistributionPointer(const int k) { //Get a pointer to the the distribution at
                                                 //lattice point k and pointing in direction 0
        return &mv_Distribution[k * stencil::Q];

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

        return &mv_OldDistribution[k * stencil::Q];

    }
    inline double* getDistributionOldPointer(const int k) {

        return &mv_OldDistribution[k * stencil::Q];

    }
    inline const double& getDistributionOld(const int k) const {

        return mv_OldDistribution[k];

    }
    inline double& getDistributionOld(const int k) {

        return mv_OldDistribution[k];

    }

    int m_Q = stencil::Q;

};

template<template<class> class obj, class lattice, typename T, int num=1> //obj template will guarantee a unique instance of the class with its own
                                       //static vector. I pass the class to itself to guarantee this
                                       //
                                       //T determines the type stored in the vector, num
                                       //determines the number of directions for each lattice point (eg you might
                                       //want two directions corresponding to velocity in x and velocity in y at
                                       //each point.
class Parameter {

    public:

        Parameter() {

            mv_Parameter.resize(num * lattice::m_N); //Resize to the desired size

        }
        Parameter(const Parameter& other) {}

        inline std::vector<T>& getParameter() const { //Returns const vector containing the parameter

            return mv_Parameter;

        }
        inline std::vector<T>& getParameter() { //Returns vector containing the parameter

            return mv_Parameter;

        }
        inline T* getParameterPointer(const int k) const {  //Returns const pointer to parameter at lattice point k and
                                                    //direction 0
            return &mv_Parameter[k * m_Num];

        }
        inline T* getParameterPointer(const int k) { //Returns pointer to parameter at lattice point k and
                                             //direction 0
            return &mv_Parameter[k * m_Num];

        }
        inline T& getParameter(const int idx) const { //Returns const parameter at index idx

            return mv_Parameter[idx];

        }
        inline T& getParameter(const int idx) { //Returns const parameter at index idx

            return mv_Parameter[idx];
            
        }

        static inline void initialise(const T val,const int k, const int idx=0){

            if (!mm_Initialised.count(k*m_Num+idx)) mv_Parameter[k*m_Num+idx]=val;
            else mm_Initialised.erase(k*m_Num+idx);

        }

        static inline void initialiseUser(const T val,const int k, const int idx=0){

            if (!mm_Initialised.count(k*m_Num+idx)) {
                mv_Parameter[k*m_Num+idx]=val;
                mm_Initialised.insert({k*m_Num+idx,true});
            }
            else mm_Initialised.erase(k*m_Num+idx);            

        }

        static std::map<int,bool> mm_Initialised;

        static constexpr int m_Num=num;
        using ParamType = T;

        inline void Save(std::string filename, int t, std::string datadir);

        static std::vector<T> mv_Parameter; //Static vector (Does not change between objects of the class)
        
};

template<template<class> class obj, class lattice, typename T, int num>
std::vector<T> Parameter<obj, lattice, T, num>::mv_Parameter; //Must allocate memory for static vector outside of class

template<template<class> class obj, class lattice, typename T, int num>
std::map<int,bool> Parameter<obj, lattice, T, num>::mm_Initialised;

template<template<class> class obj,class lattice,  typename T, int num>
inline void Parameter<obj, lattice, T, num>::Save(std::string filename, int t, std::string datadir) { //Function to save parameter stored in this class

    char fdump[512];
    sprintf(fdump, (datadir + filename + "_t%li.mat").c_str(), t); //Buffer containing file name and location.

#ifdef MPIPARALLEL //When MPI is used we need a different approach for saving as all nodes are trying to write to the file

    MPI_File fh;

    MPI_File_open(MPI_COMM_SELF, fdump, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh); //Open the file using mpi in write only mode
    
    MPI_File_seek(fh, sizeof(double) * CURPROCESSOR * m_Num * (lattice::m_LX * lattice::m_LY * lattice::m_LZ) / NUMPROCESSORS, MPI_SEEK_SET); //Skip to a certain location in the file, currently
    MPI_File_write(fh,&mv_Parameter[lattice::m_HaloSize * m_Num], m_Num * (lattice::m_N - 2 * lattice::m_HaloSize), MPI_DOUBLE, MPI_STATUSES_IGNORE);

    MPI_File_close(&fh);
         
#else

    std::ofstream fs(fdump, std::ios::out | std::ios::binary);
    
    fs.seekp(sizeof(double) * CURPROCESSOR * m_Num * (lattice::m_LX * lattice::m_LY * lattice::m_LZ) / NUMPROCESSORS);

    for (int k = lattice::m_HaloSize; k < lattice::m_N - lattice::m_HaloSize; k++) { 

        for(int idx = 0; idx < m_Num; idx++) fs.write((char *)(&mv_Parameter[k * m_Num + idx]), sizeof(double));
        
    };

    fs.close();

#endif

}

template<class lattice, template<class> class ...parameters>
class ParameterSave {

    public:

        ParameterSave(std::string datadir, int saveinterval = 1) : m_SaveInterval(saveinterval), m_DataDir(datadir) {

            int status = system(((std::string)"mkdir -p " + m_DataDir).c_str());
            if (status) std::cout << "Error creating output directory" << std::endl;

        }
        ParameterSave(ParameterSave<lattice,parameters...>& other) : m_SaveInterval(other.saveinterval), m_DataDir(other.datadir) {

            int status = system(((std::string)"mkdir -p " + m_DataDir).c_str());
            if (status) std::cout << "Error creating output directory" << std::endl;

        }
        inline void Save(int timestep);

    private:

        const int m_SaveInterval;
        std::string m_DataDir;
        std::tuple<parameters<lattice>...> mt_Parameters;
        
};

template<class lattice, template<class> class ...parameters>
inline void ParameterSave<lattice, parameters...>::Save(int timestep) {

    std::string dir = m_DataDir;

    if (timestep % m_SaveInterval == 0) {

        if (CURPROCESSOR == 0) std::cout << "SAVING at timestep " << timestep << std::endl;

        if constexpr (sizeof...(parameters) != 0) {

            std::apply([timestep,dir](parameters<lattice>&... params) {

                (params.Save(params.m_Name,timestep,dir),...);

            }, mt_Parameters);

        }
        
    }
}

template<class lattice>
struct Velocity : public Parameter<Velocity, lattice, double, lattice::m_NDIM> {

    static constexpr char m_Name[] = "Velocity";

}; //Velocity, with directions D corresponding to the number of cartesian directions in the stencilUw

template<class lattice>
struct Density : public Parameter<Density, lattice, double> {

    static constexpr char m_Name[] = "Density";

}; //Density

template<class lattice>
struct Pressure : public Parameter<Pressure, lattice, double>{

    static constexpr char m_Name[] = "Pressure";

}; //Presure

template<class lattice>
struct OrderParameter : public Parameter<OrderParameter, lattice, double> {

    static constexpr char m_Name[] = "OrderParameter";

}; //Order parameter representing relative concentration of the phases

template<class lattice>
struct ChemicalPotential : public Parameter<ChemicalPotential, lattice, double> {

    static constexpr char m_Name[] = "ChemicalPotential";

}; //Chemical potential for the multicomponent model

template<class lattice>
struct LaplacianOrderParameter : public Parameter<LaplacianOrderParameter, lattice, double> {

    static constexpr char m_Name[] = "LaplacianOrderParameter";

}; //Laplacian of the order parameter

template<class lattice>
struct GradientOrderParameter : public Parameter<GradientOrderParameter, lattice, double, lattice::m_NDIM> {

    static constexpr char m_Name[]="GradientOrderParameter";

}; //Directional first order gradients of the order parameter

template<class lattice>
struct SolidLabels : public Parameter<SolidLabels, lattice, int> {
    static constexpr char m_Name[] = "SolidLabels";
}; //Labelling of geometry
