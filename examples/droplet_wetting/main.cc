#include <lbm.hh>

// This script simulates a droplet on a flat surface with a given contact angle.


const int lx = 60; // Size of domain in x direction
const int ly = 60; // Size of domain in y direction
const int lz = 1; // Size of domain in z direction

const int timesteps = 1000000; // Number of iterations to perform
const int saveInterval = 50000; // Interval to save global data

const double contactAngle = 120; // Contact angle of the liquid on the solid
const double dropRadius = 20; // Radius to initialise the droplet


// Set up the lattice, including the resolution and data/parallelisation method
using Lattice = LatticeProperties<DataOldNew, X_Parallel<1>, lx, ly, lz>;



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
    double y0 = 1.5;// - dropRadius * cos(contactAngle*M_PI/180);
    double r2 = pow(x-lx/2.0, 2) + pow(y-y0, 2);
    if (r2 < pow(dropRadius,2)) return 1;
    else return -1;
}



// Modify the traits of the binary model to use MRT
using TraitFlowField = DefaultTraitFlowField<Lattice> ::SetCollisionOperator<MRT>;


int main(int argc, char **argv){
    mpi.init();

    // Define the models to be used
    FlowFieldBinary<Lattice,TraitFlowField> flowFieldModel; //Flowfield (navier stokes solver) that can be used with the binary model
    Binary<Lattice> componentSeparationModel; //Binary model with hybrid equilibrium and forcing term

    componentSeparationModel.getPreProcessor<ChemicalPotentialCalculatorBinary>().setA(0.015);
    componentSeparationModel.getPreProcessor<ChemicalPotentialCalculatorBinary>().setKappa(0.03);

    // Define the contact angle on the solid
    componentSeparationModel.getPreProcessor<CubicWetting>().setThetaDegrees(contactAngle);
    componentSeparationModel.getPreProcessor<CubicWetting>().setAlpha(sqrt(2));
    componentSeparationModel.setTau2(0.51);

    // Define the solid using the function above
    SolidLabels<>::set<Lattice>(initSolid);

    // Initialise the liquid and gas using the function above
    OrderParameter<>::set<Lattice>(initFluid);

    // Algorithm creates an object that combines the lattice models and runs them in order
    Algorithm lbm(flowFieldModel, componentSeparationModel);

    // Set up the handler object for saving data
    ParameterSave<Lattice> saver("data/");
    saver.SaveHeader(timesteps, saveInterval); // Create a header with lattice information (lx, ly, lz, NDIM (2D or 3D), timesteps, saveInterval)

    // Perform the main LBM loop
    for (int timestep=0; timestep<=timesteps; timestep++) {
        if (timestep%saveInterval==0) {
            std::cout<<"Saving at timestep "<<timestep<<"."<<std::endl;
            saver.SaveParameter<SolidLabels<>>(timestep);
            saver.SaveParameter<OrderParameter<>>(timestep);
            saver.SaveParameter<Velocity<>,Lattice::NDIM>(timestep);
        }
        lbm.evolve();
    }

    return 0;
}
