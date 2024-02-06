# LBM library

This library is a work in progress, so its functionality is subject to change. There may also be undiscovered bugs, so please let us know if you find anything or have any suggestions.

This is a header-only library, so to use it simply link to the 'src/' directory when compiling and include the 'lbm.hh' file within your script.
To run in parallel using MPI, ensure that the `-DMPIPARALLEL` flag is used during compilation. To run in parallel using OpenMP, ensure that the `-fopenmp` flag is used during compilation. The code supports hybrid parallelisation so these flags are not mutually exclusive.

You must use C++17 or later.

See the 'examples/' directory for scripts showing the library in use.

## Library overview

This library makes heavy use of templates to specify methods to use and to pass information about the lattice (stored within the `LatticeProperties` class).

### Models
The library is centered around various models which contain the functions to solve a given lattice Boltzmann equation (e.g. the collision, equilibrium distributions, etc.).
These are stored in the 'src/LBModels/' directory.

| Name | Description | Reference |
|:---:|:-------|:----|
|`FlowField`      | Standard LBM. Solves Navier--Stokes equations. | |
|`FlowFieldPressure`      | Standard LBM with pressure as the zeroth moment. Solves Navier--Stokes equations. | |
|`Binary`         | Two-component model. Solves Cahn--Hilliard equation. | doi:10.1007/978-3-319-44649-3 Section 9.2.2 |
|`BinaryFlowField`| Two-component model. Solves Navier--Stokes equations. | doi:10.1007/978-3-319-44649-3 Section 9.2.2 |
|`PressureLee`      | Modified two-component LBM with pressure as the zeroth moment. Extra forcing terms are included to improve stability at high density ratios and reduce spurious currents. | doi:10.1016/j.jcp.2010.07.007 |
|`BinaryLee`| Modified two-component model. Extra forcing terms are included to improve stability at high density ratios and reduce spurious currents. Solves Cahn--Hilliard equation. | doi:10.1016/j.jcp.2010.07.007 |
|`EvaporationHumidity`| Solves an advection-diffusion equation for the evolution of vapour within a gas phase. | doi:10.1103/PhysRevE.103.053307 |
|`PressureLeeHumidity`      | Modified PressureLee model with different pressure calculation for use with EvaporationHumidity | doi:10.1016/j.jcp.2010.07.007 doi:10.1103/PhysRevE.103.053307 |
|`BinaryLeeHumidity`| Modified BinaryLee for use with EvaporationHumidity | doi:10.1016/j.jcp.2010.07.007 doi:10.1103/PhysRevE.103.053307 |
|`PressureTernaryLee`      | Modified three-component LBM with pressure as the zeroth moment. Extra forcing terms are included to improve stability at high density ratios and reduce spurious currents. | doi:10.1103/PhysRevE.97.033312 |
|`BinaryTernaryLee`| Modified three-component model. Extra forcing terms are included to improve stability at high density ratios and reduce spurious currents. Solves Cahn--Hilliard equation. | doi:10.1103/PhysRevE.97.033312 |
|`PressureTernaryLeeHumidity`      | Modified PressureTernaryLee model with different pressure calculation for use with EvaporationHumidity | doi:10.1103/PhysRevE.97.033312 doi:10.1103/PhysRevE.103.053307 |
|`WellBalanced`      | Two-phase model with low spurious velocity. Solves Navier--Stokes equations. | doi:10.1063/5.0041446 |

Each model is given a traits template parameter that contains the stencil, boundary methods, collision method, any `PreProcessors` and `PostProcessors` such as gradient caclulation, and any number of forces/source terms.
The models each have a default trait, e.g. `DefaultTraitFlowField` for `FlowField`, but these can be modified if desired.

### Boundary Conditions
| Name | Description | Reference |
|:---:|:-------|:----|
|`BounceBack`| No-slip solid (halfway) | |
|`Bouzidi`| No-slip solid (arbitrary solid distance) | doi:10.1063/1.1399290 |
|`Bouzidi2`| Higher order no-slip solid (variable solid distance) | doi:10.1063/1.1399290 |
|`Convective`| Outflow condition assuming macroscopic variables are convcted out of the domain | doi:10.1103/PhysRevE.87.063301 |
|`Dirichlet`| Constant 0th moment on the wall (halfway) | doi:10.1007/978-3-319-44649-3 Section 8.5.2.1 |
|`DirichletVariable`| Constant 0th moment on the wall, can vary along the boundary (halfway) | doi:10.1007/978-3-319-44649-3 Section 8.5.2.1 |
|`ExtrapolationOutflow`| 0 gradient boundary condition | doi:10.1103/PhysRevE.87.063301 |
|`FreeSlip`| Slip solid/mirror boundary | doi:10.1007/978-3-319-44649-3 Section 5.3.7 |
|`InterpolatedDirichlet`| Constant 0th moment on the wall (variable solid distance) | doi:10.1016/j.jcp.2012.11.027 |
|`Neumann`| Constant gradient boundary condition | doi:10.1007/978-3-319-44649-3 Section 8.5.3.1 |
|`PressureOutflow`| Constant gradient boundary condition specific to pressure in the humidity model | doi:10.1103/PhysRevE.88.013304 |
|`Refill`| Refill condition for moving boundaries | doi:10.1103/PhysRevE.103.053307 |
|`VelocityInflow`| Constant 1st moment on the wall (halfway) | doi:10.1007/978-3-319-44649-3 Section 5.3.5.1 |
|`VelocityInflowVariable`| Constant 1th moment on the wall, can vary along the boundary (halfway) | doi:10.1007/978-3-319-44649-3 Section 5.3.5.1 |
|`ZouHe`| Constant density or velocity | doi:10.1063/1.869307 |

### Parameters
Values that vary across the lattice such as velocity and density are stored as `Parameter` objects.
Several functions are provided to set these values and to read them out.
They can also be passed to methods in the `SaveHandler` class in order to write them to a file during the simulation.
A list of the various parameters can be found in `src/Parameters.hh`.

## Running the Examples

If you navigate to to the `examples` folder, you will see two example folders. 
The `droplet_wetting` example involves a droplet on a solid substrate, that will reach a contact angle which you prescribe in the `main.cc` file.
The `poiseuille_flow` example involves a single component pushed by a body force between two solid plates.
The `layered_poiseuille_flow` example involves two component layers parallel to each other between two solid plates. These are pushed by a body force. You can set the relative viscosity of each component in the `main.cc` file.

Each `main.cc` file will include comments explaining what is happening.

To run these examples, navigate to the relevant folder and enter `make`. Then you can enter `./run.exe` to run the example. If you want to run using OpenMP, enter `export OMP_NUM_THREADS=[number]` but replace `[number]` with the number of threads you want. If you want to run using MPI, enter `mpirun -np [number] run.exe`.

You can then run the `plot.py` file in the terminal, which will produce some analysis plots showing the results.
