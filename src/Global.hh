#pragma once
#ifdef MPIPARALLEL
#include <mpi.h>
#endif
#include <functional>

/**
 * \file Global.hh
 * \brief This contains global parameters for the code. Eventually this should be mostly phased out.
 */

int NUMPROCESSORS = 1; //!<Number of processors (1 for a serial job, updated otherwise).
int CURPROCESSOR = 0; //!<Id of the current process (0 for a serial job).
char* MPIBUFFER; //!<Pointer to the MPI buffer.
int MPIBUFFERSIZE; //!<Size of the MPI buffer.

enum{x = 0, y = 1, z = 2};