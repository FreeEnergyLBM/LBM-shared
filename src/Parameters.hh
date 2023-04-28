#ifndef PARAMETER_HEADER
#define PARAMETER_HEADER
#include <vector>
#include<string>
#include<map>
#include<any>
#include "Global.hh"
#include <fstream>
#include <iostream>
using namespace std;

//Parameters.hh: This file details how macroscopic quantities are stored and interacted with.
//The Distribution class contains some vectors and functions to return values from these. This will be inherited
//from in the Data classes (see Data.hh). The Parameter class contains a static vector and functions to return
//the parameters stored. The reason for this static vector is to ensure that, if I use "Density" in multiple
//classes, the value will be consistent between classes. Note that, for every template configuration, a new
//static variable is created, so I pass each class to itself as a template to ensure the parameters are unique.

template<class  stencil> //Distribution must know information about the stencil as this determines the size
                         //of the vectors and the data layout
struct Distribution_Base{ //Distribution base class
    Distribution_Base(std::vector<int>& neighbors):mv_DistNeighbors(neighbors){

    }
    vector<double> mv_Distribution; //Vector that will store the distribution information
    vector<double> mv_OldDistribution; //Possibly unused vector storing the old distributions (at t_old=t_new-1)
                                       //This is required in the collision step
    std::vector<int>& mv_DistNeighbors; //Reference to vector containing neighbor information
    vector<double>& getDistribution(){ //Get a vector containing the distributions
        return mv_Distribution;
    }
    const double* getDistributionPointer(const int k) const{ //Get a constant pointer to the the distribution at
                                                             //lattice point k and pointing in direction 0
        return &mv_Distribution[k*stencil::Q];
    }
    double* getDistributionPointer(const int k){ //Get a pointer to the the distribution at
                                                 //lattice point k and pointing in direction 0
        return &mv_Distribution[k*stencil::Q];
    }
    const double& getDistribution(const int idx) const{ //Get const distribution value at a given index
        return mv_Distribution[idx];
    }
    double& getDistribution(const int idx){ //Get distribution value at a given index
        return mv_Distribution[idx];
    }
    vector<double>& getDistributionOld(){ //SAME AS ABOVE BUT FOR OLD DISTRIBUTION
        return mv_OldDistribution;
    }
    const double* getDistributionOldPointer(const int k) const{
        return &mv_OldDistribution[k*stencil::Q];
    }
    double* getDistributionOldPointer(const int k){
        return &mv_OldDistribution[k*stencil::Q];
    }
    const double& getDistributionOld(const int k) const{
        return mv_OldDistribution[k];
    }
    double& getDistributionOld(const int k){
        return mv_OldDistribution[k];
    }

    int m_Q=stencil::Q;

};

template<class obj,typename T,int num> //obj template will guarantee a unique instance of the class with its own
                                       //static vector. I pass the class to itself to guarantee this
                                       //
                                       //T determines the type stored in the vector, num
                                       //determines the number of directions for each lattice point (eg you might
                                       //want two directions corresponding to velocity in x and velocity in y at
                                       //each point.
class Parameter{
    public:

        Parameter(){
            mv_Parameter.resize(num*N); //Resize to the desired size
        }

        vector<T>& getParameter() const{ //Returns const vector containing the parameter
            return mv_Parameter;
        }
        vector<T>& getParameter(){ //Returns vector containing the parameter
            return mv_Parameter;
        }
        T* getParameterPointer(const int k) const{  //Returns const pointer to parameter at lattice point k and
                                                    //direction 0
            return &mv_Parameter[k*num];
        }
        T* getParameterPointer(const int k){ //Returns pointer to parameter at lattice point k and
                                             //direction 0
            return &mv_Parameter[k*num];
        }
        T& getParameter(const int idx) const{ //Returns const parameter at index idx
            return mv_Parameter[idx];
        }
        T& getParameter(const int idx){ //Returns const parameter at index idx
            return mv_Parameter[idx];
        }

        const static int m_Num=num;
        using ParamType=T;

        static void Save(std::string filename,int t,std::string datadir);

    private:
        static vector<T> mv_Parameter; //Static vector (Does not change between objects of the class)
        
};

template<class obj,typename T,int num>
vector<T> Parameter<obj,T,num>::mv_Parameter; //Must allocate memory for static vector outside of class

template<class obj,typename T,int num>
void Parameter<obj,T,num>::Save(std::string filename,int t,std::string datadir){ //Function to save parameter stored in this class

    char fdump[512];
    sprintf(fdump, (datadir+filename+"_t%li.mat").c_str(),t); //Buffer containing file name and location.

#ifdef MPIPARALLEL //When MPI is used we need a different approach for saving as all nodes are trying to write to the file

    MPI_File fh;

    MPI_File_open(MPI_COMM_SELF, fdump,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh); //Open the file using mpi in write only mode
    
    MPI_File_seek(fh,sizeof(double)*CURPROCESSOR*num*(LX*LY*LZ)/NUMPROCESSORS,MPI_SEEK_SET); //Skip to a certain location in the file, currently
    MPI_File_write(fh,&mv_Parameter[MAXNEIGHBORS*LY*LZ*num],num*(N-2*MAXNEIGHBORS*LY*LZ),MPI_DOUBLE,&status);
    //MPI_File_write_all(fh,&mv_Parameter[MAXNEIGHBORS*LY*LZ*num],num*(N-2*MAXNEIGHBORS*LY*LZ),MPI_DOUBLE,&status);
    //for (int k = MAXNEIGHBORS*LY*LZ; k < N-MAXNEIGHBORS*LY*LZ; k++ ) { 
    //    
    //    for(int idx=0;idx<num;idx++) {
    //        MPI_File_write(fh,&mv_Parameter[k*num+idx],num*,MPI_DOUBLE,&status);
    //    }

    //};

    MPI_File_close(&fh);
         
#else

    std::ofstream fs(fdump, std::ios::out | std::ios::binary );
    
    fs.seekp(sizeof(double)*CURPROCESSOR*num*(LX*LY*LZ)/NUMPROCESSORS);

    for (int k = MAXNEIGHBORS*LY*LZ; k < N-MAXNEIGHBORS*LY*LZ; k++ ) { 

        for(int idx=0;idx<num;idx++) fs.write((char *)(&mv_Parameter[k*num+idx]), sizeof(double));
        
    };

    fs.close();

#endif

}

template<class ...parameters>
class ParameterSave{
    public:
        ParameterSave(std::string datadir,int saveinterval=1):m_SaveInterval(saveinterval),m_DataDir(datadir){
            system(((std::string)"mkdir "+m_DataDir).c_str());
        }
        void Save(int timestep){
            if (timestep%SAVEINTERVAL==0) {
                if(CURPROCESSOR==0) std::cout<<"SAVING at timestep "<<timestep<<""<<std::endl;
                (parameters::Save(parameters::m_Name,timestep,m_DataDir),...);
            }
        }
    private:
        const int m_SaveInterval;
        const std::string m_DataDir;
};

template<typename T,int ndim>
struct Velocity : public Parameter<Velocity<T,ndim>,T,ndim>{static constexpr char m_Name[]="Velocity";}; //Velocity, with directions D
                                                                        //corresponding to the number of cartesian
                                                                        //directions in the stencilUw

template<typename T>
struct Density : public Parameter<Density<T>,T,1>{static constexpr char m_Name[]="Density";}; //Density

template<typename T>
struct Pressure : public Parameter<Pressure<T>,T,1>{static constexpr char m_Name[]="Pressure";}; //Presure

template<typename T>
struct OrderParameter : public Parameter<OrderParameter<T>,T,1>{static constexpr char m_Name[]="OrderParameter";}; //Order parameter representing relative
                                                                   //concentration of the phases

template<typename T>
struct ChemicalPotential : public Parameter<ChemicalPotential<T>,T,1>{static constexpr char m_Name[]="ChemicalPotential";}; //Chemical potential for the multicomponent model

template<typename T>
struct LaplacianOrderParameter : public Parameter<LaplacianOrderParameter<T>,T,1>{static constexpr char m_Name[]="LaplacianOrderParameter";}; //Laplacian of the order parameter

template<typename T,int ndim>
struct GradientOrderParameter : public Parameter<GradientOrderParameter<T,ndim>,T,ndim>{static constexpr char m_Name[]="GradientOrderParameter";}; //Directional first order gradients of the order parameter

struct SolidLabels : public Parameter<SolidLabels,int,1>{static constexpr char m_Name[]="SolidLabels";}; //Labelling of geometry



#endif