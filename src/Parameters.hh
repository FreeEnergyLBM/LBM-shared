#ifndef PARAMETER_HEADER
#define PARAMETER_HEADER
#include <vector>
#include<string>
#include<map>
#include<any>
#include "Global.hh"
using namespace std;
template<class  stencil,class model>
struct Distribution{

    static vector<double> mv_Distribution;
    static vector<double> mv_OldDistribution;

    vector<double>& getDistribution() const{
        return mv_Distribution;
    }
    vector<double>& getDistribution(){
        return mv_Distribution;
    }
    const double* getDistributionPointer(const int k) const{
        return &mv_Distribution[k*stencil::Q];
    }
    double* getDistributionPointer(const int k){
        return &mv_Distribution[k*stencil::Q];
    }
    const double& getDistribution(const int k) const{
        return mv_Distribution[k];
    }
    double& getDistribution(const int k){
        return mv_Distribution[k];
    }
    vector<double>& getDistributionOld() const{
        return mv_OldDistribution;
    }
    vector<double>& getDistributionOld(){
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
    virtual int streamIndex(int k,int Q){
        return 0;
    };

};

template<class  stencil,class model>
vector<double> Distribution<stencil, model>::mv_Distribution;

template<class  stencil,class model>
vector<double> Distribution<stencil, model>::mv_OldDistribution;

template<class obj,typename T,int num>
class Parameter{
    public:

        Parameter(){
            mv_Parameter.resize(num*N);
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

template<typename T>
struct Pressure : public Parameter<Pressure<T>,T,1>{};

template<typename T>
struct OrderParameter : public Parameter<OrderParameter<T>,T,1>{};

struct SolidLabels : public Parameter<SolidLabels,int,1>{};

#endif