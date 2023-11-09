#include <lbm.hh>
#include <math.h>
#include <stdlib.h>

// This script simulates a four component layered poiseuille flow setup.
// You can modify the densities and relaxation times via the 'setDensities' and 'setTaus' functions. 
// You can modify the body force magnitude in the setMagnitudeX function

int lx = 200; // Size of domain in x direction
int ly = 100; // Size of domain in y direction
int timesteps = 5000; // Number of iterations to perform
int saveInterval = 1000; // Interval to save global data
double radius = 25.0;
double theta = 90.0;
double A = 0.025;
double kappa = 0.05;
double Hsat = 0.3;
double Hwall = 0.0;
double offsetx = 0.0;
double offsety = -ly/2.+12;
double dens1 = 1;
double dens2 = 1;
double diffusivity = 0.02;
std::string datadir = "data/";

InputParameters params;

//using Lattice = LatticeProperties<ParallelX<3>, lx, ly>;
using Lattice = LatticePropertiesRuntime<ParallelX<3>, 2>;
//using Lattice = LatticeProperties<NoParallel, lx, ly>;

// Function used to define the solid geometry
// Here we set a solid at the top and bottom, in the conditions that return 1;
int initBoundary(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx-1)/2. - offsetx) * (xx - (lx-1)/2. - offsetx) + (yy - (ly-1)/2. - offsety) * (yy - (ly-1)/2. - offsety);
    //return 0;
    //if (yy > 1 && (yy >= ly - 1 || xx == 1 || xx == lx - 2)) return 4;
    if (yy > 1 && (yy >= ly - 1 || xx == 1 || xx == lx - 2)) return 4;
    if (yy <= 1) return 1;
    //if (yy == 1 || yy == ly - 2 || xx == 1 || xx == lx - 2) return 4;
    //if (xx <= 1) return 1;
    //else if (xx == lx - 2) return 4;
    //else if (xx >= lx - 1) return 1;
    if(sqrt(rr2)<radius) return 5;
    //if(sqrt(rr2)<radius) return 5;
    //else if (xx < lx/2.+offset) return 5;
    return 0;
}

double initFluid(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx-1)/2. - offsetx) * (xx - (lx-1)/2. - offsetx) + (yy - (ly-1)/2. - offsety) * (yy - (ly-1)/2. - offsety);
    //return 0.25*((double)rand()/(double)RAND_MAX);
    
    if (yy <= 1) return 0;
    else return 0.5-0.5*tanh(2*(sqrt(rr2)-radius)/(sqrt(8*kappa/A)));
    //return 0.5-0.5*tanh(2*((xx - lx/2.-offset))/(sqrt(8*kappa/A)));
    //return 0.5-0.5*tanh(2*((yy - ly/2.))/(sqrt(8*kappa/A)));
    //if(sqrt(rr2)<radius) return 1;
    //else return 0;
}

double initHumidity(int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx-1)/2. - offsetx) * (xx - (lx-1)/2. - offsetx) + (yy - (ly-1)/2. - offsety) * (yy - (ly-1)/2. - offsety);
    //return 0.25*((double)rand()/(double)RAND_MAX);
    if(sqrt(rr2)<radius&&yy > 1) return Hsat;
    else return 0;
    //if (xx <= lx/2.) return 0.5;
    //if (xx <= 0 || xx >= lx - 1) return 0;
    //return 0.5 - 0.5*(xx - (lx/2.+offset))/(lx-lx/2.-offset-1.5);
    
    //else return 0;
    //int yy = computeY(ly, 1, k);
    //return 0.5*tanh(2*(yy-3*ly/4)/(D))-0*0.5*tanh(2*(yy-ly)/(D))+0.5-0*0.5*tanh(2*(yy)/(D));
}
/*
bool interfaceCondition(const double& val, int k){
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    //return 0.25*((double)rand()/(double)RAND_MAX);
    if(sqrt(rr2)<radius) return false;
    else return true;
    //return val<0.5;
}
*/

bool interfaceCondition(const double& val, int k){
    return val<0.5;
}

bool interfaceConditionK(int k){
    return OrderParameter<>::get<Lattice>(k)<0.5;
}

// Function used to define the fluid
// Here we set a tanh transition in the y direction. This should match the equilibrium profile.
/*
double initFluid1(int k) {
    //int xx = computeXGlobal<Lattice>(k);
    //double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    //return 0.25*((double)rand()/(double)RAND_MAX);
    //return 0.5+0.5*tanh((sqrt(rr2)-radius)/(sqrt(2*kappa/A)));
    int yy = computeY(ly, 1, k);
    return 0.5*tanh(2*(yy-3*ly/4)/(D))-0*0.5*tanh(2*(yy-ly)/(D))+0.5-0*0.5*tanh(2*(yy)/(D));
}

double initFluid2(int k) {
    int yy = computeY(ly, 1, k);
    return 0.5*tanh(2*(yy-ly/2)/(D))-0.5*tanh(2*(yy-3*ly/4)/(D));
}

double initFluid3(int k) {
    int yy = computeY(ly, 1, k);
    return 0.5*tanh(2*(yy-ly/4)/(D))-0.5*tanh(2*(yy-ly/2)/(D));
}
*/

using traithumid = DefaultTraitHumidity<Lattice>::SetStencil<D2Q9>;
using traitpressure = typename DefaultTraitPressureLeeHumidity<Lattice> :: template SetBoundary<PressureOutflow<typename DefaultTraitPressureLeeHumidity<Lattice>::Forces>,BounceBack>;
//using traitpressure = typename DefaultTraitPressureLee<Lattice> :: AddForce<BodyForce<>>;

double distancefunc(int k, int idx){

    using Stencil = typename traithumid::Stencil;

    double normaldist = -0.5*sqrt(8*kappa/A)*atanh((OrderParameter<>::get<Lattice>(k)-0.5)/0.5);
    double* gradorderparam = GradientOrderParameter<>::getAddress<Lattice,Lattice::NDIM>(k,0);
    //float_or_double magnitudegradient = sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2)+pow(gradorderparam[2], 2));
    std::vector<double> normal;
    if constexpr (Lattice::NDIM==1) normal = {-gradorderparam[0] / sqrt(pow(gradorderparam[0], 2))};
    else if constexpr (Lattice::NDIM==2) normal = {-gradorderparam[0] / sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2)), -gradorderparam[1] / sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2))};
    else normal = {-gradorderparam[0] / sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2)+pow(gradorderparam[2], 2)), -gradorderparam[1] / sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2)+pow(gradorderparam[2], 2)), -gradorderparam[2] / sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2)+pow(gradorderparam[2], 2))};
    double normdotci=0;
    double magci=0;
    for (int xyz = 0; xyz<Lattice::NDIM; xyz++){
        normdotci +=  normal[xyz] * Stencil::Ci_xyz(xyz)[idx];
        magci += Stencil::Ci_xyz(xyz)[idx] * Stencil::Ci_xyz(xyz)[idx];
    }
    //return 0.5; //!!!!!!!!!!!!!!!!
    double dist = fabs(normaldist / normdotci/sqrt(magci));

    if (idx == 0 || std::isnan(dist) || std::isinf(dist)) {

        return 0.5;
    }
    else if (dist<=0) {
        return 0.5;
    }
    return dist;

}

void initParams(std::string inputfile) {
    
    params.addParameter<int>(lx,"lx");
    params.addParameter<int>(ly,"ly");
    params.addParameter<int>(timesteps,"timesteps");
    params.addParameter<int>(saveInterval,"saveInterval");
    params.addParameter<double>(radius,"radius");
    params.addParameter<double>(theta,"theta");
    params.addParameter<double>(A,"A");
    params.addParameter<double>(kappa,"kappa");
    params.addParameter<double>(Hsat,"Hsat");
    params.addParameter<double>(Hwall,"Hwall");
    params.addParameter<double>(offsetx,"offsetx");
    params.addParameter<double>(offsety,"offsety");
    params.addParameter<double>(dens1,"dens1");
    params.addParameter<double>(dens2,"dens2");
    params.addParameter<double>(diffusivity,"diffusivity");
    params.addParameter<std::string>(datadir,"datadir");

    params.readInput(inputfile);

    Lattice::init(lx,ly,1);

}

template<typename TTrait = typename traithumid::SetDataType<DataOldNewEquilibrium>>
EvaporationHumidity<Lattice, TTrait> initHumidity(){

    EvaporationHumidity<Lattice,TTrait> humidity;

    using dbtype = InterpolatedDirichlet;

    humidity.template getBoundary<InterpolatedDirichlet>().setInterfaceDistanceFunction(distancefunc);
    humidity.template getBoundary<dbtype>().setInterfaceID(5);
    humidity.template getBoundary<dbtype>().setInterfaceVal(Hsat);

    humidity.template getBoundary<Dirichlet>().setInterfaceID(4);
    humidity.template getBoundary<Dirichlet>().setInterfaceVal(Hwall);

    humidity.setDiffusivity(diffusivity);

    humidity.template getPreProcessor<HumidityBoundaryLabels>().setInterfaceCondition(interfaceCondition);
    humidity.template getPreProcessor<SetHumidityLiquid>().setInterfaceVal(Hsat);

    humidity.template getPostProcessor<Gradients<Humidity<>, CentralXYZInterfaceNoSolid>>().setInterfaceDistance(distancefunc);
    //humidity.getPostProcessor<GradientsInterface<Humidity<>, CentralXYZInterfaceNoSolid>>().setInterfaceCondition(interfaceConditionK);
    humidity.template getPostProcessor<Gradients<Humidity<>, CentralXYZInterfaceNoSolid>>().setInterfaceVal(Hsat);

    humidity.template getForce<EvaporationHumiditySource<EvaporationSourceMethod>>().setInterfaceHumidity(Hsat);
    humidity.template getForce<EvaporationHumiditySource<EvaporationSourceMethod>>().setGasDensity(dens2);

    return humidity;

}

template<typename TTrait = typename traitpressure::SetDataType<DataOldNewEquilibrium>>
PressureLeeHumidity<Lattice, TTrait> initPressure(){

    PressureLeeHumidity<Lattice,TTrait> pressure;

    pressure.template getForce<EvaporationPressureSource<EvaporationSourceMethod>>().setInterfaceHumidity(Hsat);
    pressure.template getForce<EvaporationPressureSource<EvaporationSourceMethod>>().setGasDensity(dens2);
    pressure.template getBoundary<PressureOutflow<typename DefaultTraitPressureLeeHumidity<Lattice>::Forces>>().setPressureCalculator(pressure.computePressure);
    pressure.template getBoundary<PressureOutflow<typename DefaultTraitPressureLeeHumidity<Lattice>::Forces>>().setForceTuple(pressure.mt_Forces);

    return pressure;

}

template<typename TTrait = typename DefaultTraitBinaryLeeHumidity<Lattice>:: template SetDataType<DataOldNewEquilibrium>>
BinaryLeeHumidity<Lattice, TTrait> initBinary(){

    BinaryLeeHumidity<Lattice, TTrait> binary;

    binary.setDensity1(dens1);
    binary.setDensity2(dens2);

    binary.template getForce<MuSourceLocal>().setBeta(A);
    binary.template getForce<MuSourceNonLocal>().setBeta(A);
    
    binary.template getPostProcessor<ChemicalPotentialCalculatorBinaryLee>().setA(A);
    binary.template getPostProcessor<ChemicalPotentialCalculatorBinaryLee>().setKappa(kappa);
    
    binary.template getPreProcessor<MassLossCalculatorInterpolated>().setInterfaceHumidity(Hsat);
    binary.template getPreProcessor<MassLossCalculatorInterpolated>().setDiffusivity(diffusivity);
    binary.template getPreProcessor<MassLossCalculatorInterpolated>().setInterfaceWidth(sqrt(8*kappa/A));
    binary.template getPreProcessor<MassLossCalculatorInterpolated>().setPhiGasLiquid(0,1);
    binary.template getPreProcessor<MassLossCalculatorInterpolated>().setGasDensity(dens2);
    
    double wettingprefactor = - cos(theta*M_PI/180.0)*sqrt(2*A/kappa);

    binary.template getPreProcessor<GradientsMultiStencil<OrderParameter<>, CentralXYZNoSolid, CentralQNoSolid, MixedXYZNoSolid, MixedQNoSolid, LaplacianCentralWetting>>().setWettingPrefactor(wettingprefactor);

    return binary;

}