#ifndef PARAMETER_HEADER
#define PARAMETER_HEADER
#include <vector>
#include<string>
#include<map>
#include<any>
#include "Global.hh"
using namespace std;
template<class  stencil, template<typename givenstencil> class datatype,int i>
struct Distribution{
    Distribution(){
        if constexpr (datatype<stencil>::OneArray==true){
            mv_Distribution.resize(2*stencil::Q*LX*LY*LZ);
        }
        else {
            mv_Distribution.resize(stencil::Q*LX*LY*LZ);
            mv_OldDistribution.resize(stencil::Q*LX*LY*LZ);
        }
    }
    static vector<double> mv_Distribution;
    static vector<double> mv_OldDistribution;
    vector<double>& getDistribution() const{
        return mv_Distribution;
    }
    vector<double>& getDistribution(){
        return mv_Distribution;
    }
    const double* getDistributionPointer(const int k) const{
        return &mv_Distribution[k];
    }
    double* getDistributionPointer(const int k){
        return &mv_Distribution[k];
    }
    const double& getDistribution(const int k) const{
        return mv_Distribution[k];
    }
    double& getDistribution(const int k){
        return mv_Distribution[k];
    }
    vector<double> getDistributionOld() const{
        return mv_OldDistribution;
    }
    vector<double> getDistributionOld(){
        return mv_OldDistribution;
    }
    const double* getDistributionOldPointer(const int k) const{
        return &mv_OldDistribution[k];
    }
    double* getDistributionOldPointer(const int k){
        return &mv_OldDistribution[k];
    }
    const double& getDistributionOld(const int k) const{
        return mv_OldDistribution[k];
    }
    double& getDistributionOld(const int k){
        return mv_OldDistribution[k];
    }
    
    void iterate();
    void getNeighbor();
    void communicate();
    int offset;
};

template<class  stencil, template<typename givenstencil> class datatype,int i>
vector<double> Distribution<stencil, datatype, i>::mv_Distribution;

template<class  stencil, template<typename givenstencil> class datatype,int i>
vector<double> Distribution<stencil, datatype, i>::mv_OldDistribution;

template<class obj,typename T,int num>
class Parameter{
    public:
        Parameter(){
            mv_Parameter.resize(num*LX*LY*LZ);
        }
        vector<T>& getParameter() const{
            return mv_Parameter;
        }
        vector<T>& getParameter(){
            return mv_Parameter;
        }
        T* getParameterPointer(const int k) const{
            return &mv_Parameter[k];
        }
        T* getParameterPointer(const int k){
            return &mv_Parameter[k];
        }
        T& getParameter(const int k) const{
            return mv_Parameter[k];
        }
        T& getParameter(const int k){
            return mv_Parameter[k];
        }
    private:
        static vector<T> mv_Parameter;
};

template<class obj,typename T,int num>
vector<T> Parameter<obj,T,num>::mv_Parameter;

template<typename T,typename stencil>
struct Velocity : public Parameter<Velocity<T,stencil>,T,stencil::D>{};

template<typename T>
struct Density : public Parameter<Density<T>,T,1>{};

struct SolidLabels : public Parameter<SolidLabels,int,1>{};

#endif