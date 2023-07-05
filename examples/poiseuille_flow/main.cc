#include <lbm.hh>

// This script simulates a Poiseuille flow in a channel driven by an external force


const int lx = 100; // Size of domain in x direction
const int ly = 50; // Size of domain in y direction

const int timesteps = 10000; // Number of iterations to perform
const int saveInterval = 1000; // Interval to save global data

const double force = 1e-6; // Driving force, equivalent to the pressure gradient


// Function used to define the solid geometry
// Here we set a solid at the top and bottom
int initSolid(const int k) {
    int y = computeY(ly, 1, k);
    if (y <= 1 || y >= ly - 2) return 1;
    else return 0;
}


int main(int argc, char **argv){
    mpi.init();

    // Set up the lattice, including the resolution and data/parallelisation method
    using Lattice = LatticeProperties<Data1, X_Parallel<1>, lx, ly>;

    // We use the 'FlowField' LBM model, which is the standard Navier-Stokes solver
    // We need to modify the traits of the model to include a body force as an 'AddOn'.
    using Force = BodyForce<Lattice,Guo<Lattice,D2Q9>>;
    using PoiseuilleTrait = DefaultTraitFlowField<Lattice> ::AddAddOn<Force>;
    FlowField<Lattice,PoiseuilleTrait> model;

    // Define the magnitude of the body force
    model.getAddOn<Force>().setMagnitudeX(force);

    // Define the solid using the function above
    SolidLabels<Lattice> solid;
    solid.set(initSolid);

    // Set up the handler object for saving data
    ParameterSave<Lattice,SolidLabels,Velocity> saver("data/"); // The templates specify the parameters you want to save (and the lattice)
    saver.SaveHeader(timesteps, saveInterval); // Create a header with lattice information (lx, ly, lz, NDIM (2D or 3D), timesteps, saveInterval)

    // Algorithm creates an object that can run our chosen LBM model
    Algorithm lbm(model);

    // Perform necessary initialisation for the models
    lbm.initialise();

    // Perform the main LBM loop
    for (int timestep=0; timestep<=timesteps; timestep++) {
        if (timestep%saveInterval==0) saver.Save(timestep);
        lbm.evolve();
    }

    return 0;
}
