#include <lbm.hh>
#include <math.h>
#include <stdlib.h>

// This script simulates a four component layered poiseuille flow setup.
// You can modify the densities and relaxation times via the 'setDensities' and 'setTaus' functions. 
// You can modify the body force magnitude in the setMagnitudeX function

const int lx = 10; // Size of domain in x direction
const int ly = 200; // Size of domain in y direction

const int timesteps = 40000; // Number of iterations to perform
const int saveInterval = 1000; // Interval to save global data

//Parameters to control the surface tension and width of the diffuse interface
//Use these if you want the surface tensions to all be the same
//double BETA=0.001;
//double GAMMA=-BETA*6.25;

double D = 5;
double SURFACETENSION = 0.001;

double MOBILITY = 0.01;

using Lattice = LatticeProperties<DataOldNew, NoParallel, lx, ly>;
const int NUM_COMPONENTS=4; //Number of fluid components

template<int compid>
using traitNCOMPChemPotCalculator 
    = typename DefaultTraitNComponent<Lattice, compid, NUM_COMPONENTS> 
                :: template AddPreProcessor<ChemicalPotentialCalculatorNComponent,
                                  GradientsMultiStencil< OrderParameter<NUM_COMPONENTS - 1>, CentralXYZNoSolid,LaplacianCentralNoSolid>>;

using traitNCOMPPressure = DefaultTraitFlowFieldPressureNComp<Lattice,NUM_COMPONENTS> 
                            :: template AddForce<BodyForce<>>
                            :: template SetCollisionOperator<MRT>;

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

int main(int argc, char **argv){
    //mpi.init();

    // Set up the lattice, including the resolution and data/parallelisation method
    
    // We need to modify the traits of the navier stokes model to include a bodyforce and change the collision model to MRT, so that high viscosity contrasts produce accurate results.
    // We need to add a chemical potential calculator to one of the models too, given to the first NComponent model.

    FlowFieldPressureNComp<Lattice,NUM_COMPONENTS,traitNCOMPPressure> PressureNavierStokes;
    NComponent<Lattice, 0, NUM_COMPONENTS,traitNCOMPChemPotCalculator<0>> NCompAllenCahn1;
    NComponent<Lattice, 1, NUM_COMPONENTS> NCompAllenCahn2;
    NComponent<Lattice, 2, NUM_COMPONENTS> NCompAllenCahn3;

    PressureNavierStokes.getForce<BodyForce<>>().setMagnitudeX(0.0000001);
    PressureNavierStokes.setDensities(1.0,1.0,1.0,1.0);

    std::vector<double> taus = {1.0,1.0,1.0,1.0};

    PressureNavierStokes.setTaus(taus);

    NCompAllenCahn1.setTauAndMobility(1.0,MOBILITY);
    NCompAllenCahn2.setTauAndMobility(1.0,MOBILITY);
    NCompAllenCahn3.setTauAndMobility(1.0,MOBILITY);
    
    //Set the beta and gamma parameters in the module used to calculate the chemical potential
    auto& ChemPotCalculator = NCompAllenCahn1.getPreProcessor<ChemicalPotentialCalculatorNComponent>();

    ChemPotCalculator.setD(D); // Component 0 and 1
    
    ChemPotCalculator.setSurfaceTension(0,1,SURFACETENSION); // Component 0 and 1
    ChemPotCalculator.setSurfaceTension(0,2,SURFACETENSION); // Component 0 and 2
    ChemPotCalculator.setSurfaceTension(0,3,SURFACETENSION); // etc
    ChemPotCalculator.setSurfaceTension(1,2,SURFACETENSION);
    ChemPotCalculator.setSurfaceTension(1,3,SURFACETENSION);
    ChemPotCalculator.setSurfaceTension(2,3,SURFACETENSION);

    //Feel free to change this so you can modify the interface width and sigma parameters instead

    //Set the interface width
    NCompAllenCahn1.getForce<AllenCahnSource<AllenCahnSourceMethod,0>>().setDTauAndMobility(D, 1.0, MOBILITY);
    NCompAllenCahn2.getForce<AllenCahnSource<AllenCahnSourceMethod,1>>().setDTauAndMobility(D, 1.0, MOBILITY);
    NCompAllenCahn3.getForce<AllenCahnSource<AllenCahnSourceMethod,2>>().setDTauAndMobility(D, 1.0, MOBILITY);

    // Define the solid and fluid using the functions above
    SolidLabels<>::set<Lattice>(initSolid);
    OrderParameter<NUM_COMPONENTS-1>::set<Lattice,1,0>(initFluid1);
    OrderParameter<NUM_COMPONENTS-1>::set<Lattice,1,1>(initFluid2);
    OrderParameter<NUM_COMPONENTS-1>::set<Lattice,1,2>(initFluid3);

    // Algorithm creates an object that can run our chosen LBM model
    //Algorithm lbm(PressureNavierStokes,NCompAllenCahn1,NCompAllenCahn2,NCompAllenCahn3);
    Algorithm lbm(PressureNavierStokes,NCompAllenCahn1,NCompAllenCahn2,NCompAllenCahn3);
    // Set up the handler object for saving data
    ParameterSave<Lattice> saver("data/");
    saver.SaveHeader(timesteps, saveInterval); // Create a header with lattice information (lx, ly, lz, NDIM (2D or 3D), timesteps, saveInterval)

    // Perform the main LBM loop
    for (int timestep=0; timestep<=timesteps; timestep++) {

        // Save the desired parameters, producing a binary file for each.
        if (timestep%saveInterval==0) {
            if(mpi.rank==0)std::cout<<"Saving at timestep "<<timestep<<"."<<std::endl;
            saver.SaveParameter<SolidLabels<>>(timestep);
            saver.SaveParameter<OrderParameter<NUM_COMPONENTS-1>>(timestep);
            saver.SaveParameter<Density<>>(timestep);
            saver.SaveParameter<Velocity<>,Lattice::NDIM>(timestep);
        }
        
        // Evolve by one timestep
        lbm.evolve();
        
    }
    return 0;
}
