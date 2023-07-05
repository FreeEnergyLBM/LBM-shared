#include <lbm.hh>

// This script simulates a droplet on a flat surface with a given contact angle


const int lx = 100; // Size of domain in x direction
const int ly = 100; // Size of domain in y direction
const int lz = 1; // Size of domain in z direction

const int timesteps = 100000; // Number of iterations to perform
const int saveInterval = 10000; // Interval to save global data

const double contactAngle = 135; // Contact angle of the liquid on the solid
const double dropRadius = 20; // Radius to initialise the droplet


// Set up the lattice, including the resolution and data/parallelisation method
using Lattice = LatticeProperties<Data1, X_Parallel<1>, lx, ly, lz>;



// Function used to define the solid geometry
int initSolid(const int k) {
    int y = computeY(ly, lz, k);
    if (y <= 1 || y >= ly - 2) return 1;
    else return 0;
}



// Function used to initialise the liquid (1) and gas (-1)
double initFluid(int k) {
    int x = computeXGlobal<Lattice>(k); // global function used because the x direction is split among the processors
    int y = computeY(ly, lz, k);

    // Check if within droplet radius
    double y0 = 1.5 - dropRadius * cos(contactAngle*M_PI/180);
    double r2 = pow(x-lx/2.0, 2) + pow(y-y0, 2);
    if (r2 < pow(dropRadius,2)) return 1;
    else return -1;
}



// Modify the traits of the binary model
using TraitBinary = DefaultTraitBinary<Lattice> ::SetCollisionModel<SRT>;


int main(int argc, char **argv){
    mpi.init();

    // Define the models to be used
    FlowFieldBinary<Lattice> model1; //Flowfield (navier stokes solver) that can be used with the binary model (there are nuances with this model)
    Binary<Lattice,TraitBinary> model2; //Binary model with hybrid equilibrium and forcing term

    model1.getAddOn<ChemicalForceBinary<Lattice,Guo<Lattice,typename DefaultTraitFlowFieldBinary<Lattice>::Stencil>>>().setA(0.00015);
    model1.getAddOn<ChemicalForceBinary<Lattice,Guo<Lattice,typename DefaultTraitFlowFieldBinary<Lattice>::Stencil>>>().setKappa(0.0003);

    // Define the contact angle on the solid
    model2.getAddOn<CubicWetting<Lattice,typename TraitBinary::Stencil>>().setThetaDegrees(contactAngle);

    // Define the solid using the function above
    SolidLabels<Lattice> solid;
    solid.set(initSolid);

    // Initialise the liquid and gas using the function above
    OrderParameter<Lattice> orderParam;
    orderParam.set(initFluid);

    // Algorithm creates an object that combines the lattice models and runs them in order
    Algorithm lbm(model1, model2);

    // Set up the handler object for saving data
    ParameterSave<Lattice,SolidLabels,OrderParameter,Velocity> saver("data/"); // The templates specify the parameters you want to save (and the lattice)
    saver.SaveHeader(timesteps, saveInterval); // Create a header with lattice information (lx, ly, lz, NDIM (2D or 3D), timesteps, saveInterval)

    // Perform necessary initialisation for the models in LBM
    lbm.initialise();

    // Perform the main LBM loop
    for (int timestep=0; timestep<=timesteps; timestep++) {
        if (timestep%saveInterval==0) saver.Save(timestep);
        lbm.evolve();
    }

    return 0;
}
