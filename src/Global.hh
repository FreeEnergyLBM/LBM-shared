#pragma once

/**
 * \file Global.hh
 * \brief This contains global parameters for the code. Currently it just holds parameters containing the MPI buffer and its size.
 */

char* MPIBUFFER; //!<Pointer to the MPI buffer.
int MPIBUFFERSIZE; //!<Size of the MPI buffer.

enum{x = 0, y = 1, z = 2};
