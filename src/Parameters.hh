#include <vector>
#include<string>
#include<map>
#include<any>
using namespace std;
template<typename T,template<typename giventype> class  stencil, template<typename givenstencil> class datatype>
struct Distribution{
    Distribution(){
        if constexpr (datatype<stencil<T>>::OneArray==true){
            mv_Distribution.reserve(2*stencil<T>::m_Q*LX*LY*LZ);
        }
        else {
            mv_Distribution.reserve(stencil<T>::m_Q*LX*LY*LZ);
            mv_OldDistribution.reserve(stencil<T>::m_Q*LX*LY*LZ);
        }
    }
    vector<T> mv_Distribution;
    vector<T> mv_OldDistribution;
    void iterate();
    void getNeighbor();
    void communicate();
    int offset;
    const static stencil<T> m_Stencil;
};

template<class obj,typename T,int num>
struct Parameter{
    Parameter(){
        mv_Parameter.reserve(num*LX*LY*LZ);
    }
    static vector<T> mv_Parameter;
};

template<typename T,typename stencil>
struct Velocity : public Parameter<Velocity<T,stencil>,T,stencil::m_D>{};

template<typename T,int numbercomponents>
struct OrderParameter : public Parameter<OrderParameter<T,numbercomponents>,T,numbercomponents>{};

template<typename T,typename stencil,int numbercomponents>
struct Grad_OrderParameter_xyz : public Parameter<Grad_OrderParameter_xyz<T,stencil,numbercomponents>,T,stencil::m_Q*numbercomponents>{};

template<typename T,typename stencil,int numbercomponents>
struct Grad_OrderParameter_Q : public Parameter<Grad_OrderParameter_Q<T,stencil,numbercomponents>,T,stencil::m_Q*numbercomponents>{};

template<typename T>
struct Pressure : public Parameter<Pressure<T>,T,1>{};

struct SolidLabels : public Parameter<SolidLabels,int,1>{};

template<typename T,int numbercomponents>
struct Laplacian_OrderParameter : public Parameter<Laplacian_OrderParameter<T,numbercomponents>,T,numbercomponents>{};

template<typename T,int numbercomponents>
struct ChemicalPotential : public Parameter<ChemicalPotential<T,numbercomponents>,T,numbercomponents>{};

template<typename T>
struct simulation_parameter{
    simulation_parameter(const T initial,const std::string nm):name{nm}{
        parameter=initial;
        //VarTypes[name]=typeid(T).name();
        VarKeys[name]=&parameter;
    }
    simulation_parameter(const std::string nm):name{nm}{
        //VarTypes[name]=typeid(T).name();
        VarKeys[name]=&parameter;
    }
    static std::map<std::string, T*> VarKeys;
    //static std::map<std::string, std::string> VarTypes;
    std::string name;
    T parameter;
};



float_or_double& DensityC1;
float_or_double& DensityC2;



simulation_parameter<long> nbIter_init("nbIter");
long& nbIter=nbIter_init.parameter;

simulation_parameter<int> LX_init(1,"LX");
int& LX=LX_init.parameter;

simulation_parameter<int> LY_init(1,"LY");
int& LY=LY_init.parameter;

simulation_parameter<int> LZ_init(1,"LZ");
int& LZ=LZ_init.parameter;

simulation_parameter<int> dropletR_init(0,"dropletR");
int& dropletR=dropletR_init.parameter;

simulation_parameter<int> dropletCenterX_init(0,"dropletCenterX");
int& dropletCenterX=dropletCenterX_init.parameter;

simulation_parameter<int> dropletCenterY_init(0,"dropletCenterY");
int& dropletCenterY=dropletCenterY_init.parameter;

simulation_parameter<int> dropletCenterZ_init(0,"dropletCenterZ");
int& dropletCenterZ=dropletCenterZ_init.parameter;

simulation_parameter<float_or_double> gx_init(0,"gx");
float_or_double& gx=gx_init.parameter;
simulation_parameter<float_or_double> gy_init(0,"gy");
float_or_double& gy=gy_init.parameter;
simulation_parameter<float_or_double> gz_init(0,"gz");
float_or_double& gz=gz_init.parameter;
simulation_parameter<float_or_double> gc1_init(0,"gc1");
float_or_double& gc1=gc1_init.parameter;
simulation_parameter<float_or_double> gc2_init(0,"gc2");
float_or_double& gc2=gc2_init.parameter;
simulation_parameter<float_or_double> gc3_init(0,"gc3");
float_or_double& gc3=gc3_init.parameter;

simulation_parameter<float_or_double> initUX_init(0,"initUX");
float_or_double& initUX=initUX_init.parameter;

simulation_parameter<float_or_double> initUY_init(0,"initUY");
float_or_double& initUY=initUY_init.parameter;

simulation_parameter<float_or_double> initUZ_init(0,"initUZ");
float_or_double& initUZ=initUZ_init.parameter;

simulation_parameter<float_or_double> ushear0_init(0,"ushear0");
float_or_double& ushear0=ushear0_init.parameter;

simulation_parameter<float_or_double> tau1_init(1.0,"tau1");
float_or_double& tau1=tau1_init.parameter;

simulation_parameter<float_or_double> tau2_init(1.0,"tau2");
float_or_double& tau2=tau2_init.parameter;

simulation_parameter<float_or_double> taun_init(1.0,"taun");
float_or_double& taun=taun_init.parameter;

simulation_parameter<float_or_double> taup_init(1.0,"taup");
float_or_double& taup=taup_init.parameter;

simulation_parameter<float_or_double> alpha_init(1.0,"alpha");
float_or_double& alpha=alpha_init.parameter;

simulation_parameter<float_or_double> gamma12_init("gamma12");
float_or_double& gamma12=gamma12_init.parameter;

simulation_parameter<float_or_double> theta_init(-90,"theta");
float_or_double& theta=theta_init.parameter;

simulation_parameter<long> global_info_Step_init("global_info_Step");
long& global_info_Step=global_info_Step_init.parameter;

simulation_parameter<std::string> data_dir_init("mathematica/","data_dir");
std::string& data_dir=data_dir_init.parameter;