#pragma once
#pragma once
#include <vector>
#include<string>
#include<map>
#include<any>
#include "Global.hh"
#include "Lattice.hh"
#include <fstream>
#include <iostream>

//Parameters.hh: This file details how macroscopic quantities are stored and interacted with.
//The Distribution class contains some vectors and functions to return values from these. This will be inherited
//from in the Data classes (see Data.hh). The Parameter class contains a static vector and functions to return
//the parameters stored. The reason for this static vector is to ensure that, if I use "Density" in multiple
//classes, the value will be consistent between classes. Note that, for every template configuration, a new
//static variable is created, so I pass each class to itself as a template to ensure the parameters are unique.

template<class  stencil> //Distribution must know information about the stencil as this determines the size
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

    int getOpposite(int idx) {

        return ma_Opposites[idx];

    }

    std::vector<double> mv_Distribution; //Vector that will store the distribution information
    std::vector<double> mv_OldDistribution; //Possibly unused vector storing the old distributions (at t_old=t_new-1)
                                       //This is required in the collision step
    std::vector<int>& mv_DistNeighbors; //Reference to vector containing neighbor information
    std::vector<double>& getDistribution() { //Get a vector containing the distributions

        return mv_Distribution;

    }
    const double* getDistributionPointer(const int k) const { //Get a constant pointer to the the distribution at
                                                             //lattice point k and pointing in direction 0
        return &mv_Distribution[k * stencil::Q];
    }
    double* getDistributionPointer(const int k) { //Get a pointer to the the distribution at
                                                 //lattice point k and pointing in direction 0
        return &mv_Distribution[k * stencil::Q];

    }
    const double& getDistribution(const int idx) const { //Get const distribution value at a given index

        return mv_Distribution[idx];

    }
    double& getDistribution(const int idx) { //Get distribution value at a given index

        return mv_Distribution[idx];

    }
    std::vector<double>& getDistributionOld() { //SAME AS ABOVE BUT FOR OLD DISTRIBUTION

        return mv_OldDistribution;

    }
    const double* getDistributionOldPointer(const int k) const {

        return &mv_OldDistribution[k * stencil::Q];

    }
    double* getDistributionOldPointer(const int k) {

        return &mv_OldDistribution[k * stencil::Q];

    }
    const double& getDistributionOld(const int k) const {

        return mv_OldDistribution[k];

    }
    double& getDistributionOld(const int k) {

        return mv_OldDistribution[k];

    }

    int m_Q = stencil::Q;

};

template<class obj, typename T> //obj template will guarantee a unique instance of the class with its own
                                       //static vector. I pass the class to itself to guarantee this
                                       //
                                       //T determines the type stored in the vector, num
                                       //determines the number of directions for each lattice point (eg you might
                                       //want two directions corresponding to velocity in x and velocity in y at
                                       //each point.
class Parameter {

    public:

        Parameter(int num = 1) : m_Num(num) {

            mv_Parameter.resize(m_Num * GETPROPERTIES().m_N); //Resize to the desired size

        }
        Parameter(const Parameter& other) : m_Num(other.m_Num) {

            mv_Parameter.resize(m_Num * GETPROPERTIES().m_N); //Resize to the desired size

        }
        std::vector<T>& getParameter() const { //Returns const vector containing the parameter

            return mv_Parameter;

        }
        std::vector<T>& getParameter() { //Returns vector containing the parameter

            return mv_Parameter;

        }
        T* getParameterPointer(const int k) const {  //Returns const pointer to parameter at lattice point k and
                                                    //direction 0
            return &mv_Parameter[k * m_Num];

        }
        T* getParameterPointer(const int k) { //Returns pointer to parameter at lattice point k and
                                             //direction 0
            return &mv_Parameter[k * m_Num];

        }
        T& getParameter(const int idx) const { //Returns const parameter at index idx

            return mv_Parameter[idx];

        }
        T& getParameter(const int idx) { //Returns const parameter at index idx

            return mv_Parameter[idx];
            
        }

        const int m_Num;
        using ParamType = T;

        void Save(std::string filename, int t, std::string datadir);

        static std::vector<T> mv_Parameter; //Static vector (Does not change between objects of the class)
        
};

template<class obj, typename T>
std::vector<T> Parameter<obj, T>::mv_Parameter; //Must allocate memory for static vector outside of class

template<class obj, typename T>
void Parameter<obj, T>::Save(std::string filename, int t, std::string datadir) { //Function to save parameter stored in this class

    char fdump[512];
    sprintf(fdump, (datadir + filename + "_t%li.mat").c_str(), t); //Buffer containing file name and location.

#ifdef MPIPARALLEL //When MPI is used we need a different approach for saving as all nodes are trying to write to the file

    MPI_File fh;

    MPI_File_open(MPI_COMM_SELF, fdump, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh); //Open the file using mpi in write only mode
    
    MPI_File_seek(fh, sizeof(double) * CURPROCESSOR * m_Num * (GETPROPERTIES().m_LX * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ) / NUMPROCESSORS, MPI_SEEK_SET); //Skip to a certain location in the file, currently
    MPI_File_write(fh,&mv_Parameter[GETPROPERTIES().m_HaloSize * m_Num], m_Num * (GETPROPERTIES().m_N - 2 * GETPROPERTIES().m_HaloSize), MPI_DOUBLE, MPI_STATUSES_IGNORE);

    MPI_File_close(&fh);
         
#else

    std::ofstream fs(fdump, std::ios::out | std::ios::binary);
    
    fs.seekp(sizeof(double) * CURPROCESSOR * m_Num * (GETPROPERTIES().m_LX * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ) / NUMPROCESSORS);

    for (int k = GETPROPERTIES().m_HaloSize; k < GETPROPERTIES().m_N - GETPROPERTIES().m_HaloSize; k++) { 

        for(int idx = 0; idx < m_Num; idx++) fs.write((char *)(&mv_Parameter[k * m_Num + idx]), sizeof(double));
        
    };

    fs.close();

#endif

}

template<class ...parameters>
class ParameterSave {

    public:

        ParameterSave(std::string datadir, int saveinterval = 1) : m_SaveInterval(saveinterval), m_DataDir(datadir) {

            int status = system(((std::string)"mkdir -p " + m_DataDir).c_str());
            if (status) std::cout << "Error creating output directory" << std::endl;

        }
        ParameterSave(ParameterSave<parameters...>& other) : m_SaveInterval(other.saveinterval), m_DataDir(other.datadir) {

            int status = system(((std::string)"mkdir -p " + m_DataDir).c_str());
            if (status) std::cout << "Error creating output directory" << std::endl;

        }
        void Save(int timestep);

    private:

        const int m_SaveInterval;
        std::string m_DataDir;
        std::tuple<parameters...> mt_Parameters;
};

template<class ...parameters>
void ParameterSave<parameters...>::Save(int timestep) {

    std::string dir = m_DataDir;

    if (timestep % m_SaveInterval == 0) {

        if (CURPROCESSOR == 0) std::cout << "SAVING at timestep " << timestep << std::endl;

        if constexpr (sizeof...(parameters) != 0) {

            std::apply([timestep,dir](parameters&... params) {

                (params.Save(params.m_Name,timestep,dir),...);

            }, mt_Parameters);

        }
        
    }
}

template<typename placeholder = void>
struct VelocityTemplate : public Parameter<VelocityTemplate<placeholder>, double> {

    VelocityTemplate() : Parameter<VelocityTemplate,double>(GETPROPERTIES().m_NDIM) {};
    static constexpr char m_Name[] = "Velocity";

}; //Velocity, with directions D corresponding to the number of cartesian directions in the stencilUw

typedef VelocityTemplate<> Velocity;

struct Density : public Parameter<Density, double> {

    static constexpr char m_Name[] = "Density";

}; //Density

struct Pressure : public Parameter<Pressure, double>{

    static constexpr char m_Name[] = "Pressure";

}; //Presure

struct OrderParameter : public Parameter<OrderParameter,double> {

    static constexpr char m_Name[] = "OrderParameter";

}; //Order parameter representing relative concentration of the phases

struct ChemicalPotential : public Parameter<ChemicalPotential, double> {

    static constexpr char m_Name[] = "ChemicalPotential";

}; //Chemical potential for the multicomponent model

struct LaplacianOrderParameter : public Parameter<LaplacianOrderParameter, double> {

    static constexpr char m_Name[] = "LaplacianOrderParameter";

}; //Laplacian of the order parameter

template<typename placeholder = void>
struct GradientOrderParameterTemplate : public Parameter<GradientOrderParameterTemplate<placeholder>, double> {

    GradientOrderParameterTemplate() : Parameter<GradientOrderParameterTemplate,double>(GETPROPERTIES().m_NDIM) {};
    static constexpr char m_Name[]="GradientOrderParameter";

}; //Directional first order gradients of the order parameter

typedef GradientOrderParameterTemplate<> GradientOrderParameter;

struct SolidLabels : public Parameter<SolidLabels, int> {
    static constexpr char m_Name[] = "SolidLabels";
}; //Labelling of geometry
