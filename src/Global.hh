#ifndef GLOBAL_HEADER
#define GLOBAL_HEADER
#ifdef MPIPARALLEL
#include <mpi.h>
#endif

/**
 * \file Global.hh
 * \brief This contains global parameters for the code. Eventually this should be mostly phased out.
 */

const int NO_NEIGHBOR=1; //!< Number of neighbors needed (controls the size of the MPI region).
const double DT=1.0; //!< Timestep each iteration.
const double CS2=1.0/3.0; //!< Lattice speed of sound squared.
int LX; //!< Size of lattice in x direction.
int LY; //!< Size of lattice in y direction.
int LZ; //!< Size of lattice in z direction.
const int TIMESTEPS=1000; //!< Number of timesteps.
const int SAVEINTERVAL=50000; //!< Timestep interval at which to save.
int N; //!< Number of lattice points.
int LXdiv; //!< LX/NUMPROCESSORS+number of neighbors*2.
int MAXNEIGHBORS=0; //!< Max neighbors needed by the simulation.
int NUMPROCESSORS=1; //!< Number of processors (1 for a serial job, updated otherwise).
int CURPROCESSOR=0; //!< Id of the current process (0 for a serial job).
constexpr int NDIM=2; //!< Number of cartesian directions.
#ifdef MPIPARALLEL
template<typename Stencil>
using Parallel=X_Parallel<Stencil,NO_NEIGHBOR>; //!< Chosen MPI parallelisation method when MPI enabled.
#else
using Parallel=No_Parallel; //!< Default parallelisation when MPI disabled (Just serial).
#endif
template<typename Stencil>
using DataType=Data1<Stencil,Parallel<Stencil>>;//!< This will change the "DataType" implementation, which will govern the access of non-local data
#ifdef MPIPARALLEL
char* MPIBUFFER; //!< Buffer for MPI_Isend.
int MPIBUFFERSIZE; //!< Size of the MPI buffer.
MPI_Status status; //!< Used by MPI_Isend and MPI_Irecv
MPI_Request request; //!< Used by MPI_Isend and MPI_Irecv
#endif
#endif