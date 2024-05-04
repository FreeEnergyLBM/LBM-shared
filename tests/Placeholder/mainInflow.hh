#include <math.h>
#include <stdlib.h>

#include <lbm.hh>

// This script simulates a four component layered poiseuille flow setup.
// You can modify the densities and relaxation times via the 'setDensities' and 'setTaus' functions.
// You can modify the body force magnitude in the setMagnitudeX function

int lx = 200;             // Size of domain in x direction
int ly = 200;             // Size of domain in y direction
int timesteps = 5000;     // Number of iterations to perform
int saveInterval = 1000;  // Interval to save global data
double radius = 25.0;
double theta = 90.0;
double A = 0.003;
double kappa = A * 9 / 8;
double surfacetension = sqrt(2 * A * kappa) / 6;
double interfacewidth = sqrt(8 * kappa / A);
double Hsat = 0.3;
double Hwall = 0.0;
double offsetx = 0.0;
double offsety = -ly / 2. + 12;
double dens1 = 1;
double dens2 = 1;
double diffusivity = 0.008;
std::string datadir = "data/";
double inflowmomentum = 0;
int equilibriumtimesteps = 0;
int postwidth = lx / 4;

InputParameters params;

// using Lattice = LatticeProperties<ParallelX<3>, lx, ly>;
// using Lattice = LatticePropertiesRuntime<ParallelX<3>, 2>;
// using Lattice = LatticeProperties<NoParallel, lx, ly>;
using Lattice = LatticePropertiesRuntime<NoParallel, 2>;

// Function used to define the solid geometry
// Here we set a solid at the top and bottom, in the conditions that return 1;
int initBoundary(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1) / 2. - offsetx) * (xx - (lx - 1) / 2. - offsetx) +
                 (yy - (ly - 1) / 2. - offsety) * (yy - (ly - 1) / 2. - offsety);
    // return 0;
    // if(yy==ly/4-1&&xx==lx - 4) return 7;

    if (yy <= 1) return 1;  // && ( yy >= ly - 1 || xx <= 1 || xx >= lx - 2 )) return 1;
    if (yy >= ly - 4) return 2;
    // if (yy==ly/4 && xx == lx-4) return 7;

    if (yy >= ly / 3 && (xx <= 1)) return 3;
    if (yy >= ly / 3 && (xx >= lx - 4)) return 4;

    // if ((xx <= 1)) return 3;
    // if ((xx >= lx - 4)) return 7;
    // if (yy<ly/3) return 1;
    if ((xx < postwidth / 2 || xx > lx - postwidth / 2) && yy < ly / 3) return 1;
    // if (yy <= 5) return 1;
    // if (yy == 1 || yy == ly - 2 || xx == 1 || xx == lx - 2) return 4;
    // if (xx <= 1) return 1;
    // else if (xx == lx - 2) return 4;
    // else if (xx >= lx - 1) return 1;
    // if(sqrt(rr2)<radius) return 5;
    // if(sqrt(rr2)<radius) return 5;
    if (yy < ly / 3 + offsety) return 5;
    // else if (xx < lx/2.+offset) return 5;
    return 0;
}

double initFluid(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1) / 2. - offsetx) * (xx - (lx - 1) / 2. - offsetx) +
                 (yy - (ly - 1) / 2. - offsety) * (yy - (ly - 1) / 2. - offsety);
    // return 0.25*((double)rand()/(double)RAND_MAX);

    // if (yy <= 1) return 0;
    // return 0.5-0.5*tanh(2*(sqrt(rr2)-radius)/(sqrt(8*kappa/A)));
    // return 0.5-0.5*tanh(2*((xx - lx/2.-offset))/(sqrt(8*kappa/A)));
    return 0.5 - 0.5 * tanh(2 * ((yy - ly / 3. - offsety)) / (sqrt(8 * kappa / A)));
    // if(sqrt(rr2)<radius) return 1;
    // else return 0;
}

double initVelocity(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    // double rr2 = (xx - (lx-1)/2. - offsetx) * (xx - (lx-1)/2. - offsetx) + (yy - (ly-1)/2. - offsety) * (yy -
    // (ly-1)/2. - offsety); return 0.25*((double)rand()/(double)RAND_MAX);

    // if (yy <= 1) return 0;
    // return 0.5-0.5*tanh(2*(sqrt(rr2)-radius)/(sqrt(8*kappa/A)));
    // return 0.5-0.5*tanh(2*((xx - lx/2.-offset))/(sqrt(8*kappa/A)));
    // int k = computeK<Lattice>(xx,yy,0);
    // std::cout<<k<<std::endl;
    // return inflowmomentum*pow(yy-ly/4-((double)(ly-4-ly/4))/2.0,2)/pow((double)(ly-4-ly/4)/2.0,2)-inflowmomentum;
    if (yy < ly / 3) return 0;
    // return inflowmomentum*(yy-ly+4)/((double)(-ly+4+ly/3))/dens2-inflowmomentum/dens2;
    return -inflowmomentum / dens2;
}

double initVelocityY(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    // double rr2 = (xx - (lx-1)/2. - offsetx) * (xx - (lx-1)/2. - offsetx) + (yy - (ly-1)/2. - offsety) * (yy -
    // (ly-1)/2. - offsety); return 0.25*((double)rand()/(double)RAND_MAX);

    // if (yy <= 1) return 0;
    // return 0.5-0.5*tanh(2*(sqrt(rr2)-radius)/(sqrt(8*kappa/A)));
    // return 0.5-0.5*tanh(2*((xx - lx/2.-offset))/(sqrt(8*kappa/A)));
    // int k = computeK<Lattice>(xx,yy,0);
    // std::cout<<k<<std::endl;
    return 0;
}

double initHumidity(int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - (lx - 1) / 2. - offsetx) * (xx - (lx - 1) / 2. - offsetx) +
                 (yy - (ly - 1) / 2. - offsety) * (yy - (ly - 1) / 2. - offsety);
    // return 0.25*((double)rand()/(double)RAND_MAX);
    // if(sqrt(rr2)<radius&&yy > 1) return Hsat;
    if (yy < ly / 3 + offsety)
        return Hsat;
    else
        return 0;
    // if (xx <= lx/2.) return 0.5;
    // if (xx <= 0 || xx >= lx - 1) return 0;
    // return 0.5 - 0.5*(xx - (lx/2.+offset))/(lx-lx/2.-offset-1.5);

    // else return 0;
    // int yy = computeY(ly, 1, k);
    // return 0.5*tanh(2*(yy-3*ly/4)/(D))-0*0.5*tanh(2*(yy-ly)/(D))+0.5-0*0.5*tanh(2*(yy)/(D));
}

bool interfaceCondition(const double& val, int k) { return val < 0.5; }

bool interfaceConditionK(int k) { return OrderParameter<>::get<Lattice>(k) < 0.5; }

// Function used to define the fluid
// Here we set a tanh transition in the y direction. This should match the equilibrium profile.

using traithumid = DefaultTraitHumidity<Lattice>::SetStencil<D2Q9>::template SetBoundary<
    InterpolatedDirichlet, Dirichlet, ExtrapolationOutflow, Refill<Humidity<>>, FreeSlip>;
using traitpressure =
    typename DefaultTraitPressureLeeHumidity<Lattice>::template SetBoundary<PressureOutflow, BounceBack>;
// using traitpressure = typename DefaultTraitPressureLee<Lattice> :: AddForce<BodyForce<>>;

double distancefunc(int k, int idx) {
    using Stencil = typename traithumid::Stencil;

    double normaldist = -0.5 * sqrt(8 * kappa / A) * atanh((OrderParameter<>::get<Lattice>(k) - 0.5) / 0.5);
    double* gradorderparam = GradientOrderParameter<>::getAddress<Lattice, Lattice::NDIM>(k, 0);
    // float_or_double magnitudegradient = sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1],
    // 2)+pow(gradorderparam[2], 2));
    std::vector<double> normal;
    if constexpr (Lattice::NDIM == 1)
        normal = {-gradorderparam[0] / sqrt(pow(gradorderparam[0], 2))};
    else if constexpr (Lattice::NDIM == 2)
        normal = {-gradorderparam[0] / sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2)),
                  -gradorderparam[1] / sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2))};
    else
        normal = {-gradorderparam[0] /
                      sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2) + pow(gradorderparam[2], 2)),
                  -gradorderparam[1] /
                      sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2) + pow(gradorderparam[2], 2)),
                  -gradorderparam[2] /
                      sqrt(pow(gradorderparam[0], 2) + pow(gradorderparam[1], 2) + pow(gradorderparam[2], 2))};
    double normdotci = 0;
    double magci = 0;
    for (int xyz = 0; xyz < Lattice::NDIM; xyz++) {
        normdotci += normal[xyz] * Stencil::Ci_xyz(xyz)[idx];
        magci += Stencil::Ci_xyz(xyz)[idx] * Stencil::Ci_xyz(xyz)[idx];
    }
    // return 0.5; //!!!!!!!!!!!!!!!!
    double dist = fabs(normaldist / normdotci);  //*sqrt(magci));

    if (idx == 0 || std::isnan(dist) || std::isinf(dist)) {
        return 0.5;
    } else if (dist <= 0) {
        return 0.5;
    }
    return dist;
}

void initParams(std::string inputfile) {
    params.addParameter<int>(lx, "lx");
    params.addParameter<int>(ly, "ly");
    params.addParameter<int>(timesteps, "timesteps");
    params.addParameter<int>(saveInterval, "saveInterval");
    params.addParameter<double>(radius, "radius");
    params.addParameter<double>(theta, "theta");
    params.addParameter<double>(A, "A");
    params.addParameter<double>(kappa, "kappa");
    params.addParameter<double>(surfacetension, "surfacetension");
    params.addParameter<double>(interfacewidth, "interfacewidth");
    params.addParameter<double>(Hsat, "Hsat");
    params.addParameter<double>(Hwall, "Hwall");
    params.addParameter<double>(offsetx, "offsetx");
    params.addParameter<double>(offsety, "offsety");
    params.addParameter<double>(dens1, "dens1");
    params.addParameter<double>(dens2, "dens2");
    params.addParameter<double>(diffusivity, "diffusivity");
    params.addParameter<double>(inflowmomentum, "inflowmomentum");
    params.addParameter<std::string>(datadir, "datadir");
    params.addParameter<int>(equilibriumtimesteps, "equilibriumtimesteps");
    params.addParameter<int>(postwidth, "postwidth");

    params.readInput(inputfile);

    A = 12 * surfacetension / interfacewidth;
    kappa = pow(interfacewidth, 2) / 8.0 * A;

    Lattice::init(lx, ly, 1);
}

template <typename TTrait = typename traithumid::SetDataType<DataOldNewEquilibrium>>
EvaporationHumidity<Lattice, TTrait> initHumidity() {
    EvaporationHumidity<Lattice, TTrait> humidity;

    using dbtype = InterpolatedDirichlet;

    humidity.setCollideID({0, 6});

    humidity.template getBoundary<InterpolatedDirichlet>().setInterfaceDistanceFunction(distancefunc);
    humidity.template getBoundary<dbtype>().setNodeID(5);
    humidity.template getBoundary<dbtype>().setInterfaceVal(Hsat);

    humidity.template getBoundary<Dirichlet>().setNodeID({2, 3});
    humidity.template getBoundary<Dirichlet>().setInterfaceVal(Hwall);

    humidity.template getBoundary<ExtrapolationOutflow>().setNodeID({4, 7});

    humidity.setDiffusivity(diffusivity);

    humidity.template getProcessor<HumidityBoundaryLabels, 1>().setInterfaceCondition(interfaceCondition);
    humidity.template getProcessor<SetHumidityLiquid, 1>().setInterfaceVal(Hsat);

    humidity.template getProcessor<Gradients<Humidity<>, CentralXYZInterfaceMirrorSolid>>().setInterfaceDistance(
        distancefunc);
    // humidity.getPostProcessor<GradientsInterface<Humidity<>,
    // CentralXYZInterfaceBounceBack>>().setInterfaceCondition(interfaceConditionK);
    humidity.template getProcessor<Gradients<Humidity<>, CentralXYZInterfaceMirrorSolid>>().setInterfaceVal(Hsat);

    humidity.template getForce<EvaporationHumiditySource<EvaporationSourceMethod>>().setInterfaceHumidity(Hsat);
    humidity.template getForce<EvaporationHumiditySource<EvaporationSourceMethod>>().setGasDensity(dens2);

    return humidity;
}

template <class TLattice>
using DefaultTraitPressureLee2 = typename DefaultTrait<TLattice, 2>::template SetBoundary<BounceBack>;

template <class TLattice>
using DefaultTraitPressureLeeHumidityInflow = typename DefaultTraitPressureLee<TLattice>::template AddProcessor<
    Swapper<Velocity<>, VelocityOld<>, TLattice::NDIM>, SetParameterOld<Pressure<>, PressureOld<>>,
    SetParameterOld<Density<>, DensityOld<>>, SetParameterOld<OrderParameter<>, OrderParameterOld<>>,
    SetParameterOld<ChemicalPotential<>, ChemicalPotentialOld<>>>::
    template AddForce<EvaporationPressureSource<EvaporationSourceMethod>>::template SetBoundary<std::tuple<
        std::tuple<FreeSlip>,
        std::tuple<VelocityInflowVariable, Convective2<typename DefaultTraitPressureLee<Lattice>::template AddForce<
                                               EvaporationPressureSource<EvaporationSourceMethod>>::Forces>>,
        std::tuple<BounceBack>>>;  //,std::tuple<Temp0Outflow>>>;,ExtrapolationOutflow
//:: template SetBoundary<FreeSlip,BounceBack>
//:: AddBoundary<std::tuple<VelocityInflow,Convective>>
//:: AddBoundary<std::tuple<Temp0Outflow>>;

// template<typename TTrait = typename traitpressure::SetDataType<DataOldNewEquilibrium>>
template <typename TTrait =
              typename DefaultTraitPressureLeeHumidityInflow<Lattice>::template SetDataType<DataOldNewEquilibrium>>
auto initPressure() {
    PressureLeeHumidity<Lattice, TTrait> pressure;
    // PressureLee<Lattice, typename DefaultTraitPressureLee<Lattice>:: template SetDataType<DataOldNewEquilibrium>::
    // template
    // SetBoundary<std::tuple<std::tuple<FreeSlip,BounceBack>,std::tuple<VelocityInflowVariable,Convective2<typename
    // DefaultTraitPressureLee<Lattice>::Forces>>>>> pressure; PressureLee<Lattice, TTrait> pressure;

    pressure.setCollideID({0, 5, 6});

    pressure.template getForce<EvaporationPressureSource<EvaporationSourceMethod>>().setInterfaceHumidity(Hsat);
    pressure.template getForce<EvaporationPressureSource<EvaporationSourceMethod>>().setGasDensity(dens2);
    // pressure.template getBoundary<PressureOutflow<typename
    // DefaultTraitPressureLeeHumidity<Lattice>::Forces>>().setPressureCalculator(pressure.computePressure);
    // pressure.template getBoundary<PressureOutflow<typename
    // DefaultTraitPressureLeeHumidity<Lattice>::Forces>>().setForceTuple(pressure.mt_Forces);
    pressure.template getBoundary<BounceBack, 2>().setNodeID({1});
    pressure.template getBoundary<FreeSlip, 0>().setNodeID({2});
    pressure.template getBoundary<VelocityInflowVariable, 1>().setNodeID({3});
    // pressure.template getBoundary<Temp0Outflow,1>().setNodeID({7});
    std::unordered_map<int, std::vector<double>> momentumvector;
    int xx = 1;
    for (int yy = ly / 3; yy < ly - 4; yy++) {
        int k = computeK<Lattice>(xx, yy, 0);
        // std::cout<<k<<std::endl;
        // momentumvector[k] =
        // {-inflowmomentum*pow(yy-ly/4-((double)(ly-4-ly/4))/2.0,2)/pow((double)(ly-4-ly/4)/2.0,2)+inflowmomentum,0};
        momentumvector[k] = {
            -inflowmomentum * (yy - ly + 4) / ((double)(-ly + 4 + ly / 3)) / 3.0 + inflowmomentum / 3.0, 0};
        // momentumvector[k] = {inflowmomentum/3.0,0};
    }
    xx = Lattice::LXdiv - 6 - 4;
    for (int yy = ly / 3; yy < ly - 4; yy++) {
        int k = computeK<Lattice>(xx, yy, 0);
        // std::cout<<k<<std::endl;
        // momentumvector[k] =
        // {-inflowmomentum*pow(yy-ly/4-((double)(ly-4-ly/4))/2.0,2)/pow((double)(ly-4-ly/4)/2.0,2)+inflowmomentum,0};
        momentumvector[k] = {
            -inflowmomentum * (yy - ly + 4) / ((double)(-ly + 4 + ly / 3)) / 3.0 + inflowmomentum / 3.0, 0};
        // momentumvector[k] = {inflowmomentum/3.0,0};
    }
    pressure.template getBoundary<VelocityInflowVariable, 1>().setWallVelocity(momentumvector);
    // pressure.template getBoundary<VelocityInflow,1>().setWallVelocity({inflowmomentum,0});
    pressure
        .template getBoundary<Convective2<typename DefaultTraitPressureLee<Lattice>::template AddForce<
                                  EvaporationPressureSource<EvaporationSourceMethod>>::Forces>,
                              1>()
        .setNodeID({4});
    // pressure.template getBoundary<Convective2<typename DefaultTraitPressureLee<Lattice>:: template
    // AddForce<EvaporationPressureSource<EvaporationSourceMethod>>::Forces>,1>().setVelocityCalculator(pressure.computeVelocity);
    // pressure.template getBoundary<Convective2<typename DefaultTraitPressureLee<Lattice>:: template
    // AddForce<EvaporationPressureSource<EvaporationSourceMethod>>::Forces>,1>().setForceTuple(pressure.mt_Forces);
    // pressure.template getBoundary<Temp0Outflow,2>().setNodeID({7});
    // pressure.template getPreProcessor<ConstantGradientBoundary<Pressure<>>>().setNodeID(100);

    return pressure;
}

using orderparamgradients = GradientsMultiStencil<OrderParameter<>, CentralXYZMirrorSolid, CentralQMirrorSolid,
                                                  MixedXYZMirrorSolid, MixedQMirrorSolid, LaplacianCentralWetting>;
// using orderparamgradients = GradientsMultiStencil<OrderParameter<>, CentralXYZBounceBack,
// CentralQBounceBack,MixedXYZBounceBack, MixedQBounceBack,LaplacianCentralBounceBack>;

template <class TLattice>
using DefaultTraitBinaryLeeHumidityInflow = typename DefaultTraitBinaryLee<TLattice>::template SetProcessor<
    ConvectParameterBoundary<OrderParameter<>, OrderParameterOld<>>, ConvectParameterBoundary<Density<>, DensityOld<>>,
    ConvectParameterBoundary<ChemicalPotential<>, ChemicalPotentialOld<>>, orderparamgradients,
    ChemicalPotentialCalculatorBinaryLee>  // NoFluxSolid<OrderParameter<>>,NoFluxSolid<Pressure<>>,NoFluxSolid<Density<>>>
    ::template AddProcessor<
        std::tuple<Gradients<ChemicalPotential<>, LaplacianCentralMirrorSolid>,
                   GradientsMultiStencil<Pressure<>, CentralXYZBounceBack, CentralQBounceBack, MixedXYZBounceBack,
                                         MixedQBounceBack>,
                   MassLossCalculatorInterpolated>>::template AddForce<EvaporationPhaseSource<EvaporationSourceMethod>>
    //:: template SetBoundary<std::tuple<std::tuple<Dirichlet>,std::tuple<BounceBack>>>;
    //:: template SetBoundary<std::tuple<std::tuple<Convective2<typename DefaultTraitBinaryLee<Lattice>
    //:: template AddForce<EvaporationPhaseSource<EvaporationSourceMethod>>::Forces>>,std::tuple<BounceBack>>>;
    ::template SetBoundary<std::tuple<std::tuple<Convective>, std::tuple<BounceBack>>>;

// template<typename TTrait = typename DefaultTraitBinaryLeeHumidity<Lattice>:: template
// SetDataType<DataOldNewEquilibrium>>
template <typename TTrait =
              typename DefaultTraitBinaryLeeHumidityInflow<Lattice>::template SetDataType<DataOldNewEquilibrium>>
auto initBinary() {
    BinaryLeeHumidity<Lattice, TTrait> binary;
    // BinaryLee<Lattice, TTrait> binary;

    binary.setCollideID({0, 5, 6});

    binary.setDensity1(dens1);
    binary.setDensity2(dens2);

    binary.template getForce<MuSourceLocal>().setBeta(A);
    binary.template getForce<MuSourceNonLocal>().setBeta(A);
    // testt<typename std::tuple_element<0,typename TTrait::Processors>::type>();
    binary.template getProcessor<ChemicalPotentialCalculatorBinaryLee>().setA(A);
    binary.template getProcessor<ChemicalPotentialCalculatorBinaryLee>().setKappa(kappa);

    // binary.template getPreProcessor<SimpleMassLossCalculator>().setEvaporationRate(2e-5);

    // binary.template getPostProcessor<NoFluxSolid<OrderParameter<>>>().setNodeID({2,3,4});
    // binary.template getPostProcessor<NoFluxSolid<Pressure<>>>().setNodeID({2,3,4});
    binary.template getProcessor<MassLossCalculatorInterpolated, 1>().toggleCalculate(false);
    binary.template getProcessor<MassLossCalculatorInterpolated, 1>().setInterfaceHumidity(Hsat);
    binary.template getProcessor<MassLossCalculatorInterpolated, 1>().setDiffusivity(diffusivity);
    binary.template getProcessor<MassLossCalculatorInterpolated, 1>().setInterfaceWidth(sqrt(8 * kappa / A));
    binary.template getProcessor<MassLossCalculatorInterpolated, 1>().setPhiGasLiquid(0, 1);
    binary.template getProcessor<MassLossCalculatorInterpolated, 1>().setGasDensity(dens2);

    binary.template getBoundary<BounceBack, 1>().setNodeID({1, 2, 3});
    // binary.template getBoundary<Dirichlet>().setNodeID({44});
    // binary.template getBoundary<Dirichlet>().setInterfaceVal(0);
    binary.template getBoundary<Convective>().setNodeID({4});
    // binary.template getBoundary<Convective2<typename DefaultTraitBinaryLee<Lattice>:: template
    // AddForce<EvaporationPhaseSource<EvaporationSourceMethod>>::Forces>>().setVelocityCalculator(binary.computeVelocity);
    // binary.template getBoundary<Convective2<typename DefaultTraitBinaryLee<Lattice>:: template
    // AddForce<EvaporationPhaseSource<EvaporationSourceMethod>>::Forces>>().setForceTuple(binary.mt_Forces);

    double wettingprefactor = -2 * cos(theta * M_PI / 180.0) * sqrt(2 * A / kappa);

    // binary.template getPreProcessor<GradientsMultiStencil<OrderParameter<>, CentralXYZBounceBack, CentralQBounceBack,
    // MixedXYZBounceBack, MixedQBounceBack, LaplacianCentralWetting>>().setWettingPrefactor(wettingprefactor);
    binary.template getProcessor<orderparamgradients>().setWettingPrefactor(wettingprefactor);

    return binary;
}

template <typename T>
void AfterEquilibration(int eqsteps, T& model) {
    if (eqsteps == equilibriumtimesteps)
        model.template getProcessor<MassLossCalculatorInterpolated, 1>().toggleCalculate(true);
}
