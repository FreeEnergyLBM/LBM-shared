# LBM library

This library is a work in progress, so its functionality is subject to change.

This is a header-only library, so to use it simply link to the 'src/' directory when compiling and include the 'lbm.hh' file within your script.
To run in parallel using MPI, ensure that the `-DMPIPARALLEL` flag is used during compilation.
You must use C++17 or later.

See the 'examples/' directory for scripts showing the library in use.

## Library overview

This library makes heavy use of templates to specify methods to use and to pass information about the lattice (stored within the `LatticeProperties` class).

The library is centered around various models which contain the functions to solve a given lattice Boltzmann equation (e.g. the collision, equilibrium distributions, etc.).
These are stored in the 'src/LBModels/' directory.

Each model is given a traits template parameter that contains the stencil, boundary methods, collision method, and any 'AddOns' such as forces or the calculation of additional values.
The models each have a default trait, e.g. `DefaultTraitFlowField` for `FlowField`, but these can be modified if desired.

Values that vary across the lattice such as velocity and density are stored as `Parameter` objects.
Several functions are provided to set these values and to read them out.
They can also be passed as templates to the `ParameterSave` class in order to write them to a file during the simulation.
A list of the various parameters can be found in 'src/Parameters.hh'.
