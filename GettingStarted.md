# Getting Started

This is a brief user guide on getting started with our LBM code. The code functions as a library, with all the source code contained within header files. This file serves as a tutorial on how to write a simulation from scratch. You should try to work through this guide at least once to understand how the c. In most cases, however, you should try to start from one of the examples and adapt it to your particular setup.

# Making a basic program parallelised with OpenMP only
By writing #include "<insert path to directory here>/lbm.hh" in your main.cc file, you will gain access to the library.

Right now, our main.cc looks something like

    #include <lbm.hh>

    int main(int argc, char **argv) {
        return 0;
    }

## Creating the Lattice

To start, we must define a the size of the lattice we wish to use, as well as the MPI parallelisation method we are using.
For this, we use the LatticeProperties class, which accepts several template arguments.

LatticeProperties<"Class containing parallel methods","Size in X direction","Size in Y Direction,"Size in Z direction">

Template parameters are types/integral values that can be passed to classes and can be accessed at compile time.
The LatticeProperties class will use parts of the "Class containing parallel methods" and will store the size of the simulation, which can be accessed at compile time. The more parts of the code that can be accessed at compile time, the easier the compiler will be able to optimise the code.

In order to avoid writing this out multiple times, we define an alias with the "using" keyword.

```diff
#include <lbm.hh>

+ using Lattice = LatticeProperties<NoParallel,100,100>;

int main(int argc, char **argv) {
    return 0;
}
```

Note that the NoParallel class is a dummy class which does nothing, as we dont need MPI parallelisation, but OpenMP parallelisation will be enabled by default.
Also note that we have ignored the "Size in the Z direction" parameter; this is defaulted to 1 (if you only want a 2D simulation).

Defining this Lattice class is important as it allows the models to have inforation about the size of our simulation and how much memory to allocate.

## Choosing our LBM model

Now we can use this to initialise an LBM model. Let's say we just want a simple Navier-Stokes solver. This is contained within the "FlowField" model.

```diff
#include <lbm.hh>

using Lattice = LatticeProperties<NoParallel,100,100>;

int main(int argc, char **argv) {
    
+   FlowField<Lattice> NSSolver;

    return 0;

}
```

This class contains different functions to run parts of the LBM algorithm. For instance, we can call NSSolver.collide(); which will run the collision step.
We use the Algorithm class to run the whole LBM algorithm. We can call the evolve function from this class to run one LBM timestep. 

```diff
#include <lbm.hh>

using Lattice = LatticeProperties<NoParallel,100,100>;

int main(int argc, char **argv) {
    
    FlowField<Lattice> NSSolver;

+   Algorithm lbm(NSSolver);

#   //t=0
+   lbm.evolve(); 
#   //t=1
    
    return 0;

}
```

Initialisation of the model NSSolver is handled in the constructor of the lbm class.

## Initialising Macroscopic Variables

By default, velocity is zero and density is 1 everywhere.
If we wanted to set a different value of density and velocity, we could do this in several ways. Density and velocity are unique to each lattice but are global parameters. We use the set() function to initialise them.
The parameters "Density" and "Velocity" accept a template argument which is the instance number. If you want two densities on each lattice point you can use Density<0> and Density<1>. For one instance this can be left empty.
If I want to set the two instances of density, I would write `Density<0>::set<Lattice>(...);` and `Density<1>::set<Lattice>(...);`.
The "set function" accepts three template arguments. The first is the Lattice class. The second is the number of directions the parameter has. For instance, in this case our velocity is a vector with two dimensions, so we would put 2. The third is the direction index.
If I want to set the y component of velocity, I would write `Velocity<>::set<Lattice,2,1>(...);`.

```diff
#include <lbm.hh>

using Lattice = LatticeProperties<NoParallel,100,100>;

int main(int argc, char **argv) {

+   Density<>::set<Lattice>(1.0);
+   Velocity<>::set<Lattice,2>(0.001);

    FlowField<Lattice> NSSolver;

    Algorithm lbm(NSSolver);

    lbm.evolve(); 
    
    return 0;

}
```

Now we are running 1 timestep of an LBM algorithm with a density of 1.0 and a velocity of 0.001.

## Adding Forces and Boundary Conditions

Our model is currently not doing anything interesting. Let's add a body force and some boundaries to simulate Poiseuille flow.

Adding a body force is simple. Our model class accepts an extra argument which is a Trait class. This allows us to modify our model. For instance, we can add forces boundary conditions and additional processors. We can also change the velocity stencil and collision operator.
Each model has its own default values for these. For our model, the class is DefaultTraitFlowField. To add a force, we can write DefaultTraitFlowField<Lattice>::template AddForce<BodyForce<>>.
Note that BodyForce<> accepts a template argument for the forcing method you are using. By default this is Guo forcing.
Once we add this force, we can call NSSolver.getForce<BodyForce<>>.setMagnitudeX(...)

```diff
#include <lbm.hh>

using Lattice = LatticeProperties<NoParallel,100,100>;

int main(int argc, char **argv) {

    Density<>::set<Lattice>(1.0);
    Velocity<>::set<Lattice,2>(0.001);

+   using NewTrait = DefaultTraitFlowField<Lattice>::template AddForce<BodyForce<>>;

!   FlowField<Lattice,NewTrait> NSSolver;

+   NSSolver.getForce<BodyForce<>>.setMagnitudeX(0.00001);

    Algorithm lbm(NSSolver);

    lbm.evolve(); 
    
    return 0;

}
```
 
By default, the edges of the simulation domain are periodic. If we wish to include other types of boundary conditions, we need to label the parts of the simulation domain which will represent the boundary. We use the Geometry class as an interface to manipulate the labelling of the boundaries. The function Geometry<Lattice>::initialiseBoundaries(....) will do this for us. This function accepts two arguments. The first is a function which accepts an integer k (corresponding to the index of the lattice node) and then returns an integer which is the label we want on that node. An example of a possible implementation of this function is below.

    int initBoundary(const int k) {
        int yy = computeY(ly, 1, k); //Accepts arguments length in y direction, length in z direction and index
        if (yy == 0 || yy == ly - 1) return 1;
        else return 0;
    }

This sets the geometry labels to 1 on the top and bottom of nodes of the simulation, and 0 elsewhere.
The second argument of this function is a vector containing the integer labels corresponding to the fluid part of the domain. We want this to be {0}.

Now, our simulation looks like this:

```diff
#include <lbm.hh>

int lx = 100;
int ly = 100;

using Lattice = LatticeProperties<NoParallel,lx,ly>;


+ int initBoundary(const int k) {
+    int yy = computeY(ly, 1, k); //Accepts arguments length in y direction, length in z direction and index
+    if (yy == 0 || yy == ly - 1) return 1;
+    else return 0;
+ }


int main(int argc, char **argv) {

    Density<>::set<Lattice>(1.0);
    Velocity<>::set<Lattice,2>(0.001);

+   Geometry<Lattice>::initialiseBoundaries(initBoundary,{0})

    using NewTrait = DefaultTraitFlowField<Lattice>::template AddForce<BodyForce<>>;

    FlowField<Lattice,NewTrait> NSSolver;

    NSSolver.getForce<BodyForce<>>.setMagnitudeX(0.00001);

    Algorithm lbm(NSSolver);

    lbm.evolve(); 
    
    return 0;

}
```

We still need to choose the boundary conditions, so we modify our trait class with ::template SetBoundary<BounceBack>. This boundary condition enforces a velocity of 0 tangentially to the wall. Then we must set the nodes the boundary condition will apply to with NSSolver.template getBoundary<BounceBack>().setNodeID(...).

```diff
#include <lbm.hh>

int lx = 100;
int ly = 100;

using Lattice = LatticeProperties<NoParallel,lx,ly>;

int initBoundary(const int k) {
    int yy = computeY(ly, 1, k); //Accepts arguments length in y direction, length in z direction and index
    if (yy == 0 || yy == ly - 1) return 1;
    else return 0;
}

int main(int argc, char **argv) {

    Density<>::set<Lattice>(1.0);
    Velocity<>::set<Lattice,2>(0.001);

    Geometry<Lattice>::initialiseBoundaries(initBoundary,{0})

!   using NewTrait = DefaultTraitFlowField<Lattice>::template AddForce<BodyForce<>>::template SetBoundary<BounceBack>;

    FlowField<Lattice,NewTrait> NSSolver;

    NSSolver.getForce<BodyForce<>>.setMagnitudeX(0.00001);
+   NSSolver.template getBoundary<BounceBack>().setNodeID(1)

    Algorithm lbm(NSSolver);

    lbm.evolve(); 
    
    return 0;

}
```

Now we are simulating 2D Poiseuille flow over one timestep.

## Running the Simulation and Saving Data

The final step to this setup is to save the output data and add the simulation loop. We must create an instance of the SaveHandler class, which provides multiple saving functionalities. We will use the saveHeader(...) function, which saves info about the simulation size, saving interval and number of timesteps. We will also use the saveBoundaries(...) function to save our geometry labels and the saveParameter(...) function to save the density and velocity

```diff
#include <lbm.hh>

int lx = 100;
int ly = 100;

using Lattice = LatticeProperties<NoParallel,lx,ly>;

int initBoundary(const int k) {
    int yy = computeY(ly, 1, k); //Accepts arguments length in y direction, length in z direction and index
    if (yy == 0 || yy == ly - 1) return 1;
    else return 0;
}

int main(int argc, char **argv) {

    Density<>::set<Lattice>(1.0);
    Velocity<>::set<Lattice,2>(0.001);

    Geometry<Lattice>::initialiseBoundaries(initBoundary,{0})

    using NewTrait = DefaultTraitFlowField<Lattice>::template AddForce<BodyForce<>>::template SetBoundary<BounceBack>;

    FlowField<Lattice,NewTrait> NSSolver;

    NSSolver.getForce<BodyForce<>>.setMagnitudeX(0.00001);
    NSSolver.template getBoundary<BounceBack>().setNodeID(1)

    Algorithm lbm(NSSolver);

+   std::string datadirectory = "data";
+   int timesteps = 25000
+   int saveInterval = 1000

+   SaveHandler<Lattice> saver(datadirectory); //Create saving class

+   saver.saveHeader(timesteps, saveInterval); //Save simulation information

+   for (int timestep=0; timestep<=timesteps; timestep++) { //Loop over timesteps
#       // Save the desired parameters, producing a binary file for each.
+       if (timestep%saveInterval==0) {
+           if(mpi.rank==0)std::cout<<"Saving at timestep "<<timestep<<"."<<std::endl;
+           saver.saveBoundaries(timestep); //Save geometry labels
+           saver.saveParameter<Density<>>(timestep); //Save Density
+           saver.saveParameter<Velocity<>,Lattice::NDIM>(timestep); // Save Velocity. If saving a parameter with multiple directions, pass the number of directions as a second template argument.
+       }
        lbm.evolve(); //Evolve one timestep
+   }
    
    return 0;

}
```

These saving functions will create a binary file for each saved parameter.

## Visualising the Output Data

# Making a basic program parallelised with MPI
