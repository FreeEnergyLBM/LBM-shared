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

//extern std::function<void(void)> GETPROPERTIES;

//inline auto& GETPROPERTIES(); //!<Function that must be defined in your main file to return the lattice properies.

/**
 * \brief Function that will return an object of the chosen lattice properties class.
 * The class providing the lattice properties must be provided in the template parameters, then optional arguments can be specified in the function call.
 * This class contains a static member object of the lattice properties.
 * \param a List of optional arguments needed by property type.
 * \tparam type Type of lattice properties class.
 * \tparam args Types of optional arguments.
 * \return Object containing info for lattice properties.
 */
//template<class type, typename ...args>
//inline auto& getGlobal(args... a) {
//
//    static auto PROPERTIES = *new type(a...);
//    return PROPERTIES;
//    
//}