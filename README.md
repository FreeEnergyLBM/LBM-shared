# LBM library

This library is a work in progress, so its functionality is subject to change. There may also be undiscovered bugs, so please let us know if you find anything or have any suggestions.

This is a header-only library, so to use it simply link to the 'src/' directory when compiling and include the 'lbm.hh' file within your script.
To run in parallel using MPI, ensure that the `-DMPIPARALLEL` flag is used during compilation. To run in parallel using OpenMP, ensure that the `-fopenmp` flag is used during compilation. The code supports hybrid parallelisation so these flags are not mutually exclusive.

You must use C++17 or later.

See the 'examples/' directory for scripts showing the library in use.

## Library overview

This library makes heavy use of templates to specify methods to use and to pass information about the lattice (stored within the `LatticeProperties` class).

The library is centered around various models which contain the functions to solve a given lattice Boltzmann equation (e.g. the collision, equilibrium distributions, etc.).
These are stored in the 'src/LBModels/' directory.

Currently the code has a single component and a two component model (Section 9.2.2 The Lattice Boltzmann Method: Principles and Practice, T. Kruger et al. (2016)).

Each model is given a traits template parameter that contains the stencil, boundary methods, collision method, any `PreProcessors` and `PostProcessors` such as gradient caclulation, and any number of forces/source terms.
The models each have a default trait, e.g. `DefaultTraitFlowField` for `FlowField`, but these can be modified if desired.

Values that vary across the lattice such as velocity and density are stored as `Parameter` objects.
Several functions are provided to set these values and to read them out.
They can also be passed as templates to the `ParameterSave` class in order to write them to a file during the simulation.
A list of the various parameters can be found in `src/Parameters.hh`.

## Running the Examples

If you navigate to to the `examples` folder, you will see two example folders. 
The `droplet_wetting` example involves a droplet on a solid substrate, that will reach a contact angle which you prescribe in the `main.cc` file.
The `poiseuille_flow` example involves a single component pushed by a body force between two solid plates.
The `layered_poiseuille_flow` example involves two component layers parallel to each other between two solid plates. These are pushed by a body force. You can set the relative viscosity of each component in the `main.cc` file.

Each `main.cc` file will include comments explaining what is happening.

To run these examples, navigate to the relevant folder and enter `make`. Then you can enter `./run.exe` to run the example. If you want to run using OpenMP, enter `export OMP_NUM_THREADS=[number]` but replace `[number]` with the number of threads you want. If you want to run using MPI, enter `mpirun -np [number] run.exe`.

You can then run the `plot.py` file in the terminal, which will produce some analysis plots showing the results.