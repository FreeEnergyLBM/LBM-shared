#ifndef GLOBAL_HEADER
#define GLOBAL_HEADER
#ifdef MPIPARALLEL
#include <mpi.h>
#endif
/**
 * \file Global.hh
 * \brief This contains global parameters for the code. Eventually this should be mostly phased out.
 */

int MAXNEIGHBORS=0; //!< Max neighbors needed by the simulation.
int NUMPROCESSORS=1; //!< Number of processors (1 for a serial job, updated otherwise).
int CURPROCESSOR=0; //!< Id of the current process (0 for a serial job).

#include "Parallel.hh"
#include "Data.hh"
#ifdef MPIPARALLEL
template<typename Stencil>
using ParallelType=X_Parallel<Stencil,1>; //!< Chosen MPI parallelisation method when MPI enabled.
#else
using ParallelType=No_Parallel; //!< Default parallelisation when MPI disabled (Just serial).
#endif
template<typename Stencil>
using DataType=Data1<Stencil,ParallelType<Stencil>>;//!< This will change the "DataType" implementation, which will govern the access of non-local data

#endif