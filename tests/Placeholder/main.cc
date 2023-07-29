#include <lbm.hh>
#include <math.h>
#include <stdlib.h>

// This script simulates a Poiseuille flow in a channel with two immiscible layers of fluid driven by an external force. You can modify tau1 and tau2 to see how the velocity profile changes. tau1 and tau2 should be less than 10 or the accuracy degrades. They should never be less than 0.5, and values very close to 0.5 are prone to instability. 


const int lx = 10; // Size of domain in x direction
const int ly = 200; // Size of domain in y direction

const int timesteps = 50000; // Number of iterations to perform
const int saveInterval = 1000; // Interval to save global data

//Parameters to control the surface tension and width of the diffuse interface
//Surface tension (in lattice units) = sqrt(8*kappa*A/9)
//Interface width (in lattice units) = sqrt(kappa/A)
//For reference, the model used is in section 9.2.2 of The Lattice Boltzmann Method: Principles and Practice, T. Kruger et al. (2016)
double A=0.00015;
double kappa=A*3.125;//=0.0003;

const double RADIUS=25.;
using Lattice = LatticeProperties<Data1, NoParallel, lx, ly>;
const int NUM_COMPONENTS=4; //Number of fluid components

template<int compid>
using traitNCOMPChemPotCalculator 
    = typename DefaultTraitNComponent<Lattice, compid, NUM_COMPONENTS> 
                :: template AddPreProcessor<ChemicalPotentialCalculatorNComponent,
                                  GradientsMultiStencil< OrderParameter<NUM_COMPONENTS - 1>, CentralXYZ,LaplacianCentral>>;

using traitNCOMPPressure = DefaultTraitFlowFieldPressure<Lattice,NUM_COMPONENTS> 
                            :: template SetForce<PressureForceNComp<He>,BodyForce<>>
                            :: template AddPostProcessor<TauCalculatorNComp<NUM_COMPONENTS>>
                            :: template SetCollisionModel<MRT>;

// Function used to define the solid geometry
// Here we set a solid at the top and bottom, in the conditions that return 1;
int initSolid(const int k) {
    int y = computeY(ly, 1, k);
    if (y <= 1 || y >= ly - 2) return 1;
    else return 0;
}

// Function used to define the fluid
// Here we set a tanh transition in the y direction. This should match the equilibrium profile.


double initFluid1(int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    //return 0.5*((double)rand()/(double)RAND_MAX);;
    double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    //return 0.5+0.5*tanh((sqrt(rr2)-RADIUS)/(sqrt(2*kappa/A)));
    //if(rr2<=RADIUS*RADIUS) return 1;
    if(yy>=3*ly/4) return 1;
    else return 0;
}

double initFluid2(int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    //return 0.25*((double)rand()/(double)RAND_MAX);;
    double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    //return 0.5+0.5*tanh((sqrt(rr2)-RADIUS)/(sqrt(2*kappa/A)));
    //if(rr2>RADIUS*RADIUS&&yy<ly/2) return 1;
    if(yy>=ly/2&&yy<3*ly/4) return 1;
    else return 0;
}

double initFluid3(int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);
    //return 0.25*((double)rand()/(double)RAND_MAX);;
    double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    //return 0.5+0.5*tanh((sqrt(rr2)-RADIUS)/(sqrt(2*kappa/A)));
    //if(rr2>RADIUS*RADIUS&&yy>=ly/2&&xx>lx/2) return 0;
    if(yy>=ly/4&&yy<ly/2) return 1;
    else return 0;
}

int main(int argc, char **argv){
    //mpi.init();

    // Set up the lattice, including the resolution and data/parallelisation method
    

    // We need to modify the traits of the model to include a body force as an 'AddOn'.
    // We modify the default traits for the 'FlowFieldBinary' model, adding a bodyforce and setting the collision model to MRT, which improves accuracy at higher viscosity ratios

    FlowFieldPressure<Lattice,traitNCOMPPressure> PressureNavierStokes;
    NComponent<Lattice, 0, NUM_COMPONENTS,traitNCOMPChemPotCalculator<0>> NCompAllenCahn1;
    NComponent<Lattice, 1, NUM_COMPONENTS,traitNCOMPChemPotCalculator<1>> NCompAllenCahn2;
    NComponent<Lattice, 2, NUM_COMPONENTS,traitNCOMPChemPotCalculator<2>> NCompAllenCahn3;

    PressureNavierStokes.getForce<BodyForce<>>().setMagnitudeX(0.0000005);
    PressureNavierStokes.getPostProcessor<TauCalculatorNComp<NUM_COMPONENTS>>().setTaus(0.505,0.75,1.2,4.0);

    PressureNavierStokes.setTauMin(0.52);
    PressureNavierStokes.setTauMax(5.0);

    double** tempBeta = new double*[NUM_COMPONENTS];
    for(int i = 0; i < NUM_COMPONENTS; ++i)
        tempBeta[i] = new double[NUM_COMPONENTS];

    //For two component
    tempBeta[0][0]=0;
    tempBeta[1][1]=0;
    tempBeta[2][2]=0;
    tempBeta[3][3]=0;

    tempBeta[0][1]=A;
    tempBeta[0][2]=A;
    tempBeta[0][3]=A;
    tempBeta[1][0]=A;
    tempBeta[1][2]=A;
    tempBeta[1][3]=A;
    tempBeta[2][0]=A;
    tempBeta[2][1]=A;
    tempBeta[2][3]=A;
    tempBeta[3][0]=A;
    tempBeta[3][1]=A;
    tempBeta[3][2]=A;
    

    //For three component, construct a 3x3 matrix of beta (extend for more components)


    double** tempGamma = new double*[NUM_COMPONENTS];
    for(int i = 0; i < NUM_COMPONENTS; ++i)
        tempGamma[i] = new double[NUM_COMPONENTS];

    //For two component

    tempGamma[0][0]=0;
    tempGamma[1][1]=0;
    tempGamma[2][2]=0;
    tempGamma[3][3]=0;

    tempGamma[0][1]=2*kappa;
    tempGamma[0][2]=2*kappa;
    tempGamma[0][3]=2*kappa;
    tempGamma[1][0]=2*kappa;
    tempGamma[1][2]=2*kappa;
    tempGamma[1][3]=2*kappa;
    tempGamma[2][0]=2*kappa;
    tempGamma[2][1]=2*kappa;
    tempGamma[2][3]=2*kappa;
    tempGamma[3][0]=2*kappa;
    tempGamma[3][1]=2*kappa;
    tempGamma[3][2]=2*kappa;

    //For three component, construct a 3x3 matrix of beta (extend for more components)
    
    //Set the beta and gamma parameters in the module used to calculate the chemical potential (probably dont touch this)
    NCompAllenCahn1.getPreProcessor<ChemicalPotentialCalculatorNComponent>().setA(tempBeta);
    NCompAllenCahn1.getPreProcessor<ChemicalPotentialCalculatorNComponent>().setKappa(tempGamma);

    NCompAllenCahn2.getPreProcessor<ChemicalPotentialCalculatorNComponent>().setA(tempBeta);
    NCompAllenCahn2.getPreProcessor<ChemicalPotentialCalculatorNComponent>().setKappa(tempGamma);

    NCompAllenCahn3.getPreProcessor<ChemicalPotentialCalculatorNComponent>().setA(tempBeta);
    NCompAllenCahn3.getPreProcessor<ChemicalPotentialCalculatorNComponent>().setKappa(tempGamma);

    //Feel free to change this so you can modify the interface width and sigma parameters instead

    //Set the interface width
    NCompAllenCahn1.getForce<AllenCahnSource<AllenCahnSourceMethod,0>>().setAlpha(sqrt(4*2*kappa/A));
    NCompAllenCahn2.getForce<AllenCahnSource<AllenCahnSourceMethod,1>>().setAlpha(sqrt(4*2*kappa/A));
    NCompAllenCahn3.getForce<AllenCahnSource<AllenCahnSourceMethod,2>>().setAlpha(sqrt(4*2*kappa/A));

    // Define the solid and fluid using the functions above
    SolidLabels<>::set<Lattice>(initSolid);
    OrderParameter<NUM_COMPONENTS-1>::set<Lattice,1,0>(initFluid1);
    OrderParameter<NUM_COMPONENTS-1>::set<Lattice,1,1>(initFluid2);
    OrderParameter<NUM_COMPONENTS-1>::set<Lattice,1,2>(initFluid3);

    // Algorithm creates an object that can run our chosen LBM model
    Algorithm lbm(PressureNavierStokes,NCompAllenCahn1,NCompAllenCahn2,NCompAllenCahn3);

    // Set up the handler object for saving data
    ParameterSave<Lattice> saver("data/");
    saver.SaveHeader(timesteps, saveInterval); // Create a header with lattice information (lx, ly, lz, NDIM (2D or 3D), timesteps, saveInterval)

    // Perform the main LBM loop
    for (int timestep=0; timestep<=timesteps; timestep++) {

        // Save the desired parameters, producing a binary file for each.
        if (timestep%saveInterval==0) {
            std::cout<<"Saving at timestep "<<timestep<<"."<<std::endl;
            saver.SaveParameter<SolidLabels<>>(timestep);
            saver.SaveParameter<OrderParameter<NUM_COMPONENTS-1>>(timestep);
            saver.SaveParameter<ChemicalPotential<NUM_COMPONENTS>>(timestep);
            saver.SaveParameter<Velocity<>,Lattice::NDIM>(timestep);
        }
        
        // Evolve by one timestep
        lbm.evolve();
        
    }

    delete[] tempBeta;
    delete[] tempGamma;

    return 0;
}
