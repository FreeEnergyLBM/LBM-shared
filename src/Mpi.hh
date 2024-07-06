#pragma once
#ifdef MPIPARALLEL
#include <mpi.h>
#endif
#include <cassert>
#include <complex>

/**
 * \brief This class and its associated 'mpi' object provide a simple interface to access the MPI
 * rank and size and to automate MPI initialisation and finalisation.
 * No additional instances of this class should be created.
 */
class MpiClass {
   public:
    int rank = 0;
    int size = 1;

    void init() { init(nullptr, nullptr); }

    void init(int* argc, char*** argv) {
        if (initialised) return;
#ifdef MPIPARALLEL
        MPI_Init(argc, argv);
        getSizeRank(MPI_COMM_WORLD);
#endif
        initialised = true;
        initialisedHere = true;
    }

#ifdef MPIPARALLEL
    void init(MPI_Comm comm) {
        if (initialised) return;
        getSizeRank(comm);
        initialised = true;
    }
#endif

#ifdef MPIPARALLEL
    void getSizeRank(MPI_Comm comm) {
        MPI_Comm_size(comm, &size);
        MPI_Comm_rank(comm, &rank);
    }
#endif

    ~MpiClass() {
        if (!initialisedHere) return;
#ifdef MPIPARALLEL
#pragma omp master
        { MPI_Finalize(); }
#endif
    }

    MpiClass() { throw std::runtime_error("Error: Do not create new instances of MpiClass"); };
    MpiClass(int){};

    void barrier() {
#ifdef MPIPARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

   private:
    bool initialised = false;
    bool initialisedHere = false;
};

MpiClass mpi(0);
