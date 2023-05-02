#pragma once
#ifdef MPIPARALLEL
#include <mpi.h>
#endif
/**
 * \file Global.hh
 * \brief This contains global parameters for the code. Eventually this should be mostly phased out.
 */

int NUMPROCESSORS=1; //!< Number of processors (1 for a serial job, updated otherwise).
int CURPROCESSOR=0; //!< Id of the current process (0 for a serial job).
