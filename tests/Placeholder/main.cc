#include <lbm.hh>
#include <math.h>

// This script simulates a Poiseuille flow in a channel with two immiscible layers of fluid driven by an external force. You can modify tau1 and tau2 to see how the velocity profile changes. tau1 and tau2 should be less than 10 or the accuracy degrades. They should never be less than 0.5, and values very close to 0.5 are prone to instability. 


const int lx = 100; // Size of domain in x direction
const int ly = 100; // Size of domain in y direction

const int timesteps = 2; // Number of iterations to perform
const int saveInterval = 1; // Interval to save global data

//Parameters to control the surface tension and width of the diffuse interface
//Surface tension (in lattice units) = sqrt(8*kappa*A/9)
//Interface width (in lattice units) = sqrt(kappa/A)
//For reference, the model used is in section 9.2.2 of The Lattice Boltzmann Method: Principles and Practice, T. Kruger et al. (2016)
double A=0.00015;
double kappa=0.0003;

const double RADIUS=25.;
using Lattice = LatticeProperties<Data1, X_Parallel<1>, lx, ly>;
const int NUM_COMPONENTS=2; //Number of fluid components

template<int compid>
using traitNCOMPChemPotCalculator 
    = typename DefaultTraitNComponent<Lattice, compid, NUM_COMPONENTS> 
                :: template AddPreProcessor<ChemicalPotentialCalculatorNComponent,
                                  GradientsMultiStencil< OrderParameter<NUM_COMPONENTS - 1>, CentralXYZ,LaplacianCentral>>;

//using traitNCOMPPressure = DefaultTraitFlowFieldPressure<Lattice> 
//                            :: template SetForce<PressureForce<He, NUM_COMPONENTS,1>>;

// Function used to define the solid geometry
// Here we set a solid at the top and bottom, in the conditions that return 1;
int initSolid(const int k) {
    int y = computeY(ly, 1, k);
    if (y <= 1 || y >= ly - 2) return 1;
    else return 0;
}

// Function used to define the fluid
// Here we set a tanh transition in the y direction. This should match the equilibrium profile.
double initFluid(int k) {
    int xx = computeXGlobal<Lattice>(k);
    int yy = computeY(ly, 1, k);

    double rr2 = (xx - lx/2.) * (xx - lx/2.) + (yy - lx/2.) * (yy - lx/2.);
    return 0.5+0.5*tanh((sqrt(rr2)-RADIUS)/(sqrt(2*kappa/A)));
}

int main(int argc, char **argv){
    mpi.init();

    // Set up the lattice, including the resolution and data/parallelisation method
    

    // We need to modify the traits of the model to include a body force as an 'AddOn'.
    // We modify the default traits for the 'FlowFieldBinary' model, adding a bodyforce and setting the collision model to MRT, which improves accuracy at higher viscosity ratios

    //FlowFieldPressure<Lattice,traitNCOMPPressure> PressureNavierStokes;
    NComponent<Lattice, 1, NUM_COMPONENTS,traitNCOMPChemPotCalculator<1>> NCompAllenCahn1;

    double** tempBeta = new double*[NUM_COMPONENTS];
    for(int i = 0; i < NUM_COMPONENTS; ++i)
        tempBeta[i] = new double[NUM_COMPONENTS];

    //For two component
    tempBeta[0][0]=0;
    tempBeta[1][1]=0;
    tempBeta[0][1]=A;
    tempBeta[1][0]=A;

    //For three component, construct a 3x3 matrix of beta (extend for more components)


    double** tempGamma = new double*[NUM_COMPONENTS];
    for(int i = 0; i < NUM_COMPONENTS; ++i)
        tempGamma[i] = new double[NUM_COMPONENTS];

    //For two component
    tempGamma[0][0]=0;
    tempGamma[1][1]=0;
    tempGamma[0][1]=2*kappa;
    tempGamma[1][0]=2*kappa;

    //For three component, construct a 3x3 matrix of beta (extend for more components)
    
    //Set the beta and gamma parameters in the module used to calculate the chemical potential (probably dont touch this)
    NCompAllenCahn1.getPreProcessor<ChemicalPotentialCalculatorNComponent>().setA(tempBeta);
    NCompAllenCahn1.getPreProcessor<ChemicalPotentialCalculatorNComponent>().setKappa(tempGamma);
    //Feel free to change this so you can modify the interface width and sigma parameters instead

    //Set the interface width
    NCompAllenCahn1.getForce<AllenCahnSource<AllenCahnSourceMethod,1>>().setAlpha(sqrt(4*2*kappa/A));

    // Define the solid and fluid using the functions above
    SolidLabels<>::set<Lattice>(initSolid);
    OrderParameter<NUM_COMPONENTS-1>::set<Lattice>(initFluid);

    // Algorithm creates an object that can run our chosen LBM model
    Algorithm lbm(NCompAllenCahn1);

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
            saver.SaveParameter<Velocity<>,Lattice::NDIM>(timestep);
        }
        
        // Evolve by one timestep
        lbm.evolve();
        
    }

    return 0;
}
