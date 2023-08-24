#include <lbm.hh>
#include <math.h>
#include <stdlib.h>

// This script simulates a four component layered poiseuille flow setup.
// You can modify the densities and relaxation times via the 'setDensities' and 'setTaus' functions. 
// You can modify the body force magnitude in the setMagnitudeX function

const int lx = 99; // Size of domain in x direction
const int ly = 99; // Size of domain in y direction

const int timesteps = 10; // Number of iterations to perform
const int saveInterval = 1; // Interval to save global data

//Parameters to control the surface tension and width of the diffuse interface
//Use these if you want the surface tensions to all be the same
//double BETA=0.001;
//double GAMMA=-BETA*6.25;

double RADIUS = 15.0;

double A = 0.00025;
double kappa = 0.0005;

using Lattice = LatticeProperties<DataOldNewEquilibrium, NoParallel, lx, ly>;

// Function used to define the solid geometry
// Here we set a solid at the top and bottom, in the conditions that return 1;
int initBoundary(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    
    if (yy <= 0 || yy >= ly - 1 || xx <= 0 || xx >= lx - 1) return 4;
    else if(sqrt(rr2)<RADIUS) return 5;
    //else if (xx < lx/2.-0.6) return 5;
    else return 0;
}

double initFluid(const int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    //return 0.25*((double)rand()/(double)RAND_MAX);
    return 0.5-0.5*tanh(2*(sqrt(rr2)-RADIUS)/(sqrt(8*kappa/A)));
    //return 0.5-0.5*tanh(2*((xx - lx/2.+0.6))/(sqrt(8*kappa/A)));
    //if(sqrt(rr2)<RADIUS) return 1;
    //else return 0;
}

double initHumidity(int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    //return 0.25*((double)rand()/(double)RAND_MAX);
    if(sqrt(rr2)<RADIUS) return 0.5;
    //if (xx < lx/2.-0.6) return 0.8;
    else return 0;
    //int yy = computeY(ly, 1, k);
    //return 0.5*tanh(2*(yy-3*ly/4)/(D))-0*0.5*tanh(2*(yy-ly)/(D))+0.5-0*0.5*tanh(2*(yy)/(D));
}
/*
bool interfaceCondition(const double& val, int k){
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    //return 0.25*((double)rand()/(double)RAND_MAX);
    if(sqrt(rr2)<RADIUS) return false;
    else return true;
    //return val<0.5;
}
*/

bool interfaceCondition(const double& val, int k){
    return val<0.5;
}

// Function used to define the fluid
// Here we set a tanh transition in the y direction. This should match the equilibrium profile.
/*
double initFluid1(int k) {
    //int xx = computeXGlobal<Lattice>(k);
    //double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    //return 0.25*((double)rand()/(double)RAND_MAX);
    //return 0.5+0.5*tanh((sqrt(rr2)-RADIUS)/(sqrt(2*kappa/A)));
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
using traitpressure = DefaultTraitPressureLeeHumidity<Lattice> :: template AddBoundary<PressureOutflow<typename DefaultTraitPressureLeeHumidity<Lattice>::Forces>>;

double distancefunc(int idx, int k){

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
    
    double dist = fabs(normaldist * normdotci/magci/magci);
    //if(dist == 0.6) std::cout<<k<<std::endl;
    if (idx == 0 || std::isnan(dist) || std::isinf(dist)) return 0.5;
    else if (dist>=1) return 0.99;
    return dist;

}

int main(int argc, char **argv){
    //mpi.init();

    // Set up the lattice, including the resolution and data/parallelisation method
    
    // We need to modify the traits of the navier stokes model to include a bodyforce and change the collision model to MRT, so that high viscosity contrasts produce accurate results.
    // We need to add a chemical potential calculator to one of the models too, given to the first NComponent model.

    EvaporationHumidity<Lattice,traithumid> humidity;
    PressureLeeHumidity<Lattice,traitpressure> pressure;
    BinaryLeeHumidity<Lattice> binary;

   // BinaryLee<Lattice> binary;

    binary.getPostProcessor<ChemicalPotentialCalculatorBinaryLee>().setA(A);
    binary.getPostProcessor<ChemicalPotentialCalculatorBinaryLee>().setKappa(kappa);

    binary.getPreProcessor<MassLossCalculatorInterpolated>().setInterfaceHumidity(0.5);
    binary.getPreProcessor<MassLossCalculatorInterpolated>().setDiffusivity(0.2);
    binary.getPreProcessor<MassLossCalculatorInterpolated>().setInterfaceWidth(sqrt(8*kappa/A));
    binary.getPreProcessor<MassLossCalculatorInterpolated>().setPhiGasLiquid(0,1);

    using dbtype = InterpolatedDirichlet;

    //humidity.getBoundary<InterpolatedDirichlet,0>().setInterfaceDistanceFunction(distancefunc);
    //humidity.getBoundary<dbtype,0>().setInterfaceID(5);
    //humidity.getBoundary<dbtype,0>().setInterfaceVal(0.8);

    //humidity.getBoundary<InterpolatedDirichlet>().setInterfaceDistanceFunction(distancefunc);
    humidity.getBoundary<dbtype>().setInterfaceID(5);
    humidity.getBoundary<dbtype>().setInterfaceVal(0.5);

    humidity.getBoundary<Dirichlet>().setInterfaceID(1);
    humidity.getBoundary<Dirichlet>().setInterfaceVal(0.0);

    humidity.setDiffusivity(0.2);

    humidity.getPreProcessor<HumidityBoundaryLabels>().setInterfaceCondition(interfaceCondition);
    humidity.getPreProcessor<SetHumidityLiquid>().setInterfaceVal(0.5);

    humidity.getForce<EvaporationHumiditySource<EvaporationSourceMethod>>().setInterfaceHumidity(0.5);

    pressure.getForce<EvaporationPressureSource<EvaporationSourceMethod>>().setInterfaceHumidity(0.5);
    pressure.getBoundary<PressureOutflow<typename DefaultTraitPressureLeeHumidity<Lattice>::Forces>>().setPressureCalculator(pressure.computePressure);
    pressure.getBoundary<PressureOutflow<typename DefaultTraitPressureLeeHumidity<Lattice>::Forces>>().setForceTuple(pressure.mt_Forces);

    ParameterSave<Lattice> saver("data/");
    saver.SaveHeader(timesteps, saveInterval); // Create a header with lattice information (lx, ly, lz, NDIM (2D or 3D), timesteps, saveInterval)

    // Define the solid and fluid using the functions above
    Geometry<Lattice>::initialiseBoundaries(initBoundary);
    OrderParameter<>::set<Lattice>(initFluid);
    OrderParameterOld<>::set<Lattice>(initFluid);
    Humidity<>::set<Lattice>(initHumidity);

    // Algorithm creates an object that can run our chosen LBM model
    //Algorithm lbm(PressureNavierStokes,NCompAllenCahn1,NCompAllenCahn2,NCompAllenCahn3);
    
    // Set up the handler object for saving data
    
    //Algorithm lbm(humidity,pressure,binary);
    Algorithm lbm(binary,pressure,humidity);
    //Algorithm lbm(humidity);

    // Perform the main LBM loop
    for (int timestep=0; timestep<=timesteps; timestep++) {

        // Save the desired parameters, producing a binary file for each.
        if (timestep%saveInterval==0) {
            if(mpi.rank==0)std::cout<<"Saving at timestep "<<timestep<<"."<<std::endl;
            //saver.SaveParameter<BoundaryLabels<>>(timestep);
            saver.SaveBoundaries(timestep);
            saver.SaveParameter<Humidity<>>(timestep);
            saver.SaveParameter<ChemicalPotential<>>(timestep);
            saver.SaveParameter<Density<>>(timestep);
            saver.SaveParameter<Pressure<>>(timestep);
            saver.SaveParameter<OrderParameter<>>(timestep);
            saver.SaveParameter<MassSink<>>(timestep);
            saver.SaveParameter<Velocity<>,Lattice::NDIM>(timestep);
            saver.SaveParameter<VelocityOld<>,Lattice::NDIM>(timestep);
            saver.SaveParameter<GradientHumidity<>,Lattice::NDIM>(timestep);
        }
        
        // Evolve by one timestep
        lbm.evolve();
        
    }
    return 0;
}
