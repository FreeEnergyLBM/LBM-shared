#pragma once

#include <stdexcept>

#include "Global.hh"
#include "Mpi.hh"
#include "Parameters.hh"
#include "Service.hh"
/**
 * \file  Parallel.hh
 * \brief This contains classes to control the MPI parallelisation of the code.
 * The file contains a base class that will perform initialisation of the global MPI parameters. Other classes
 * inherit from this and provide functions to communicate between processes.
 */

/**
 * \brief This class simply initialises global MPI parameters when constructed.
 * MaxNeighbors is updated depending on the chosen number of neighbors. LXdiv (LX for each parallel block of
 * lattice points) is set based on the number of processors and number of neighbors chosen.
 */
template <class TDerived, int TWidth = 1>
class Parallel {
   public:
    /**
     * \brief Constructor that updates global parameters. This is contained within a class as I may move stuff
     *        from ParallelX to here.
     * MaxNeighbors is updated depending on the chosen number of neighbors. LXdiv (LX for each parallel block of
     * lattice points) is set based on the number of processors and number of neighbors chosen. TLattice::N is then
     * calculated as LXdiv*TLattice::LY*TLattice::LZ.
     */
    Parallel();

    /**
     * \brief Function to fill halos of adjacent processors with the chosen parameter adjacent to the edge.
     * \param obj Object of chosen parameter.
     */
    template <class TLattice, class TParameter>
    inline void communicateParameter(TParameter& obj);

    /**
     * \brief Function to update the vector containing the parameters in the communication region before communication
     * (implemented along X and Y). \param obj Object of chosen parameter.
     */
    template <class TLattice, class TParameter>
    inline void updateParameterBeforeCommunication(TParameter& obj);

    /**
     * \brief Function to update the vector containing the parameters in the communication region after communication
     * (implemented along X and Y). \param obj Object of chosen parameter.
     */
    template <class TLattice, class TParameter>
    inline void updateParameterAfterCommunication(TParameter& obj);

    /**
     * \brief Function to update the vector containing the distribution in the communication region (implemented along
     * X). \param obj Object of the distribution.
     */
    template <class TLattice, class TDistribution>
    inline void communicateDistribution(TDistribution& obj);

    /**
     * \brief Function to update unknown distributions in the adjacent processors streamed from the edge.
     * \param obj Object of the distribution.
     */
    template <class TLattice, class TDistribution>
    inline void updateDistributionBeforeCommunication(TDistribution& obj);

    /**
     * \brief Function to update the vector containing the distribution for the communication region (implemented along
     * X). \param obj Object of the distribution.
     */
    template <class TLattice, class TDistribution>
    inline void updateDistributionAfterCommunication(TDistribution& obj);

    template <class TLattice, class TDistribution>
    inline void updateDistributionBeforeCommunicationAll(TDistribution& obj);

    template <class TLattice, class TDistribution>
    inline void updateDistributionAfterCommunicationAll(TDistribution& obj);

    template <class TLattice, class TDistribution>
    inline void updateDistributionBeforeCommunicationAllOld(TDistribution& obj);

    template <class TLattice, class TDistribution>
    inline void updateDistributionAfterCommunicationAllOld(TDistribution& obj);

    template <class TLattice, class TDistribution>
    inline void updateDistributionBeforeCommunicationAllEquilibrium(TDistribution& obj);

    template <class TLattice, class TDistribution>
    inline void updateDistributionAfterCommunicationAllEquilibrium(TDistribution& obj);

    template <class TLattice, class TDistribution>
    inline void communicateDistributionAll(TDistribution& obj);

    template <class TLattice>
    static void destroyDistributionType() {
#ifdef MPIPARALLEL
// After you're done using the MPI datatype
#pragma omp master
        {
            if (isCustomDistributionTypeCreated && !isCustomDistributionTypeDestroyed)
                MPI_Type_free(&mDistributionType);
            isCustomDistributionTypeDestroyed = true;
        }
#pragma omp barrier
#endif
    }

    static constexpr int Width = TWidth;

   protected:
    int mMaxWidth = 0;
    // int mNumDirections = 0; //!<Number of communication channels. Must be redefined in child classes.
    std::vector<int> mNeighbors;    //!< IDs of the neighboring processes.
    std::vector<int> mI0Send;       //!< Lattice index to begin sending to each processor.
    std::vector<int> mI0Recv;       //!< Lattice index to begin receiving from each processor.
    std::vector<int> mI0SendDistr;  //!< Lattice index to begin sending to each processor for the distributions.
    std::vector<int> mI0RecvDistr;  //!< Lattice index to begin receiving from each processor for the distributions.

#ifdef MPIPARALLEL
    static MPI_Datatype mDistributionType;
    static bool isCustomDistributionTypeCreated;
    static bool isCustomDistributionTypeDestroyed;
#endif
};

#ifdef MPIPARALLEL
template <class TDerived, int TWidth>
MPI_Datatype Parallel<TDerived, TWidth>::mDistributionType = NULL;

template <class TDerived, int TWidth>
bool Parallel<TDerived, TWidth>::isCustomDistributionTypeCreated = false;

template <class TDerived, int TWidth>
bool Parallel<TDerived, TWidth>::isCustomDistributionTypeDestroyed = false;
#endif

template <class TDerived, int TWidth>
Parallel<TDerived, TWidth>::Parallel() {
    if (mMaxWidth < TWidth) mMaxWidth = TWidth;
}

/**
 * \details This will communicate the chosen parameter using MPI_Isend and MPI_Irecv, which are non-blocking methods of
 *          communication. This means that each process does not need to wait for the other processes to communicate. At
 *          the end of this function, we have a MPI_Waitall call, to ensure all processes are synced.
 */
template <class TDerived, int TWidth>
template <class TLattice, class TParameter>
inline void Parallel<TDerived, TWidth>::communicateParameter(TParameter& obj) {
#ifdef MPIPARALLEL
    if (mpi.size == 1) return;
// std::cout<<TParameter::mName<<" "<<obj.getCommParameter()[0]<<std::endl;
// exit(1);
#pragma omp master
    {
        int nNeighbors = mNeighbors.size();
        MPI_Request commrequest[2 * nNeighbors];
        MPI_Status commstatus[2 * nNeighbors];

        for (int iNeighbor = 0; iNeighbor < nNeighbors; iNeighbor++) {
            int tag = iNeighbor;
            MPI_Isend(&obj.getCommParameter()[mI0Send[iNeighbor] * obj.mNum],
                      TWidth * TLattice::Face[iNeighbor / 2] * obj.mNum,
                      mpi_get_type<typename TParameter::ParamType, TLattice>(), mNeighbors[iNeighbor], tag,
                      MPI_COMM_WORLD, &commrequest[2 * iNeighbor]);

            tag = (iNeighbor % 2 == 0) ? iNeighbor + 1 : iNeighbor - 1;
            MPI_Irecv(&obj.getCommParameter()[mI0Recv[iNeighbor] * obj.mNum],
                      TWidth * TLattice::Face[iNeighbor / 2] * obj.mNum,
                      mpi_get_type<typename TParameter::ParamType, TLattice>(), mNeighbors[iNeighbor], tag,
                      MPI_COMM_WORLD, &commrequest[2 * iNeighbor + 1]);
        }

        MPI_Waitall(2 * nNeighbors, commrequest, commstatus);
    }
#endif
}

template <class TDerived, int TWidth>
template <class TLattice, class TParameter>
inline void Parallel<TDerived, TWidth>::updateParameterBeforeCommunication(TParameter& obj) {
#pragma omp master
    {
        int lz = TLattice::LZdiv;
        int ly = TLattice::LYdiv;
        int lx = TLattice::LXdiv;
        int ln = 0;  // number of elements in the communication region already handled
        if (TLattice::HaloXWidth) {
            int lw = TLattice::HaloXWidth;
            for (int x = 0; x < 2 * lw; ++x)  // 4*lw
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < lz; ++z) {
                        int k = z + y * lz + x * ly * lz;
                        int xOffset = (x < lw) ? lw : lx - 3 * lw;  // (x<2*lw) ? 0 : lx-4*lw
                        int kGlobal = z + y * lz + (x + xOffset) * ly * lz;
                        for (int direction = 0; direction < TParameter::mDirections; direction++) {
                            obj.getCommParameter()[k * TParameter::mNum + (TParameter::mDirections > 1) * direction] =
                                obj.getParameter()[kGlobal * TParameter::mNum +
                                                   (TParameter::mDirections > 1) * direction];
                        }
                    }
            ln += 4 * lw * ly * lz;
        }
        if (TLattice::HaloYWidth) {
            int lw = TLattice::HaloYWidth;
            for (int x = 0; x < lx; ++x)  // 4*lw
                for (int y = 0; y < 2 * lw; ++y)
                    for (int z = 0; z < lz; ++z) {
                        int k = z + x * lz + y * lx * lz;
                        k += ln;
                        int yOffset = (y < lw) ? lw : ly - 3 * lw;  // (x<2*lw) ? 0 : lx-4*lw
                        int kGlobal = z + (y + yOffset) * lz + x * ly * lz;
                        for (int direction = 0; direction < TParameter::mDirections; direction++) {
                            obj.getCommParameter()[k * TParameter::mNum + (TParameter::mDirections > 1) * direction] =
                                obj.getParameter()[kGlobal * TParameter::mNum +
                                                   (TParameter::mDirections > 1) * direction];
                        }
                    }
            ln += 4 * lw * lx * lz;
        }
        if (TLattice::HaloZWidth) {
            int lw = TLattice::HaloZWidth;
            for (int x = 0; x < lx; ++x)  // 4*lw
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < 2 * lw; ++z) {
                        int k = y + x * ly + z * lx * ly;
                        k += ln;
                        int zOffset = (z < lw) ? lw : lz - 3 * lw;  // (x<2*lw) ? 0 : lx-4*lw
                        int kGlobal = (z + zOffset) + y * lz + x * ly * lz;
                        for (int direction = 0; direction < TParameter::mDirections; direction++) {
                            obj.getCommParameter()[k * TParameter::mNum + (TParameter::mDirections > 1) * direction] =
                                obj.getParameter()[kGlobal * TParameter::mNum +
                                                   (TParameter::mDirections > 1) * direction];
                        }
                    }
            ln += 4 * lw * lx * ly;
        }
    }
}

template <class TDerived, int TWidth>
template <class TLattice, class TParameter>
inline void Parallel<TDerived, TWidth>::updateParameterAfterCommunication(TParameter& obj) {
#pragma omp master
    {
        int lz = TLattice::LZdiv;
        int ly = TLattice::LYdiv;
        int lx = TLattice::LXdiv;
        int ln = 0;  // number of elements in the communication region already handled
        if (TLattice::HaloXWidth) {
            int lw = TLattice::HaloXWidth;
            for (int x = 0; x < 2 * lw; ++x)  // 4*lw
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < lz; ++z) {
                        int k = 2 * lw * ly * lz + z + y * lz + x * ly * lz;  // 2*lw*ly*lz +
                        int xOffset = (x < lw) ? 0 : lx - 2 * lw;             // x<2*lw
                        int kGlobal = z + y * lz + (x + xOffset) * ly * lz;
                        for (int direction = 0; direction < TParameter::mDirections; direction++) {
                            obj.getParameter()[kGlobal * TParameter::mNum + (TParameter::mDirections > 1) * direction] =
                                obj.getCommParameter()[k * TParameter::mNum +
                                                       (TParameter::mDirections > 1) * direction];
                        }
                    }
            ln += 4 * lw * ly * lz;
        }
        if (TLattice::HaloYWidth) {
            int lw = TLattice::HaloYWidth;
            for (int x = 0; x < lx; ++x)  // 4*lw
                for (int y = 0; y < 2 * lw; ++y)
                    for (int z = 0; z < lz; ++z) {
                        int k = 2 * lw * lx * lz + z + x * lz + y * lx * lz;  // 2*lw*ly*lz +
                        int yOffset = (y < lw) ? 0 : ly - 2 * lw;             // x<2*lw
                        int kGlobal = z + (y + yOffset) * lz + x * ly * lz;
                        for (int direction = 0; direction < TParameter::mDirections; direction++) {
                            obj.getParameter()[kGlobal * TParameter::mNum + (TParameter::mDirections > 1) * direction] =
                                obj.getCommParameter()[k * TParameter::mNum +
                                                       (TParameter::mDirections > 1) * direction];
                        }
                    }
            ln += 4 * lw * lx * lz;
        }
        if (TLattice::HaloZWidth) {
            int lw = TLattice::HaloZWidth;
            for (int x = 0; x < lx; ++x)  // 4*lw
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < 2 * lw; ++z) {
                        int k = 2 * lw * ly * lx + y + x * ly + z * lx * ly;  // 2*lw*ly*lz +
                        int zOffset = (z < lw) ? 0 : lz - 2 * lw;             // x<2*lw
                        int kGlobal = (z + zOffset) + y * lz + x * ly * lz;
                        for (int direction = 0; direction < TParameter::mDirections; direction++) {
                            obj.getParameter()[kGlobal * TParameter::mNum + (TParameter::mDirections > 1) * direction] =
                                obj.getCommParameter()[k * TParameter::mNum +
                                                       (TParameter::mDirections > 1) * direction];
                        }
                    }
            ln += 4 * lw * lx * ly;
        }
    }
}

template <class TDerived, int TWidth>
template <class TLattice, class TDistribution>
inline void Parallel<TDerived, TWidth>::communicateDistribution(TDistribution& obj) {
#ifdef MPIPARALLEL
    if (mpi.size == 1) return;

    using TStencil = typename TDistribution::Stencil;

#pragma omp master
    {
        MPI_Request commdist_request[20];

        int id = 0;

        for (int dir = 0; dir < TDerived::mNumDirections; dir++)

            for (int idx = 1; idx < TStencil::Q; idx++) {
                int direction =
                    TStencil::Ci_xyz(TDerived::mCommDirection[dir])[idx];  // 0 --> TDerived::mCommDirection[dir]
                int iNeighbor = 2 * dir;
                int tag = idx;  // iNeighbor;

                if (direction == -1) {
                    MPI_Isend(&obj.getCommDistribution()[mI0SendDistr[0] * TStencil::Q + idx], 1, mDistributionType,
                              mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commdist_request[id++]);
                    MPI_Irecv(&obj.getCommDistribution()[mI0RecvDistr[1] * TStencil::Q + idx], 1, mDistributionType,
                              mNeighbors[iNeighbor + 1], tag, MPI_COMM_WORLD, &commdist_request[id++]);

                } else if (direction == 1) {
                    MPI_Isend(&obj.getCommDistribution()[mI0SendDistr[1] * TStencil::Q + idx], 1, mDistributionType,
                              mNeighbors[iNeighbor + 1], tag, MPI_COMM_WORLD, &commdist_request[id++]);
                    MPI_Irecv(&obj.getCommDistribution()[mI0RecvDistr[0] * TStencil::Q + idx], 1, mDistributionType,
                              mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commdist_request[id++]);
                }
            }

        MPI_Waitall(id, commdist_request, MPI_STATUSES_IGNORE);
    }
#endif
}

template <class TDerived, int TWidth>
template <class TLattice, class TDistribution>
inline void Parallel<TDerived, TWidth>::communicateDistributionAll(TDistribution& obj) {
#ifdef MPIPARALLEL
    if (mpi.size == 1) return;

    using TStencil = typename TDistribution::Stencil;

#pragma omp master
    {
        int nNeighbors = mNeighbors.size();
        MPI_Request commrequest[2 * nNeighbors];
        MPI_Status commstatus[2 * nNeighbors];

        for (int iNeighbor = 0; iNeighbor < nNeighbors; iNeighbor++) {
            int tag = iNeighbor;
            MPI_Isend(&obj.getCommDistribution()[mI0SendDistr[iNeighbor] * TStencil::Q * TWidth],
                      TWidth * TLattice::LY * TLattice::LZ * TStencil::Q, mpi_get_type<double, TLattice>(),
                      mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2 * iNeighbor]);

            tag = (iNeighbor % 2 == 0) ? iNeighbor + 1 : iNeighbor - 1;
            // int iNeighbor2 = (iNeighbor % 2 == 0) ? iNeighbor + 1 : iNeighbor - 1;
            MPI_Irecv(&obj.getCommDistribution()[mI0RecvDistr[iNeighbor] * TStencil::Q * TWidth],
                      TWidth * TLattice::LY * TLattice::LZ * TStencil::Q, mpi_get_type<double, TLattice>(),
                      mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2 * iNeighbor + 1]);
        }

        MPI_Waitall(2 * nNeighbors, commrequest, commstatus);
    }
#endif
}

template <class TDerived, int TWidth>
template <class TLattice, class TDistribution>
inline void Parallel<TDerived, TWidth>::updateDistributionBeforeCommunication(TDistribution& obj) {
#pragma omp master
    {
        using TStencil = typename TDistribution::Stencil;
        int lz = TLattice::LZdiv;
        int ly = TLattice::LYdiv;
        int lx = TLattice::LXdiv;
        int lw;
        if (TLattice::HaloXWidth) {
            lw = TLattice::HaloXWidth;
            for (int x = 0; x < 2; ++x)  // x<4
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = z + y * lz + x * ly * lz;
                            int xOffset = (x < 1) ? lw - 1 : lx - lw - 1;  //(x<2) ? lw-1 : lx-lw-3
                            int kGlobal = z + y * lz + (x + xOffset) * ly * lz;
                            // std::cerr<<obj.getCommDistribution().size()<<" "<<k*TStencil::Q + idx<<std::endl;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getDistribution()[kGlobal * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloYWidth) {
            lw = TLattice::HaloYWidth;
            for (int y = 0; y < 2; ++y)  // y<4
                for (int x = 0; x < lx; ++x)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = z + x * lz + y * lx * lz;
                            int yOffset = (y < 1) ? lw - 1 : ly - lw - 1;  // (y<2) ? lw-1 : ly-lw-3
                            int kGlobal = z + (y + yOffset) * lz + x * ly * lz;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getDistribution()[kGlobal * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloZWidth) {
            lw = TLattice::HaloZWidth;
            for (int z = 0; z < 2; ++z)  // z<4
                for (int x = 0; x < lx; ++x)
                    for (int y = 0; y < ly; ++y)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = y + x * ly + z * lx * ly;
                            int zOffset = (z < 1) ? lw - 1 : lz - lw - 1;  // (z<2) ? lw-1 : lz-lw-3
                            int kGlobal = (z + zOffset) + y * lz + x * ly * lz;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getDistribution()[kGlobal * TStencil::Q + idx];
                        }
        }
    }
}

template <class TDerived, int TWidth>
template <class TLattice, class TDistribution>
inline void Parallel<TDerived, TWidth>::updateDistributionAfterCommunication(TDistribution& obj) {
#pragma omp master
    {
        using TStencil = typename TDistribution::Stencil;
        int lz = TLattice::LZdiv;
        int ly = TLattice::LYdiv;
        int lx = TLattice::LXdiv;
        int lw;
        if (TLattice::HaloXWidth) {
            lw = TLattice::HaloXWidth;
            for (int x = 0; x < 2; ++x)  // x<4
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            if ((x < 1 && TStencil::Ci_x[idx] > 0) || (x >= 1 && TStencil::Ci_x[idx] < 0)) {
                                int k = 2 * ly * lz + z + y * lz + x * ly * lz;  // 2*ly*lz +
                                int xOffset = (x < 1) ? lw : lx - lw - 2;        // (x<2) ? lw-1 : lx-lw-3
                                int kGlobal = z + y * lz + (x + xOffset) * ly * lz;
                                obj.getDistribution()[kGlobal * TStencil::Q + idx] =
                                    obj.getCommDistribution()[k * TStencil::Q + idx];
                            }
                        }
        }
        if (TLattice::HaloYWidth) {
            lw = TLattice::HaloYWidth;
            for (int y = 0; y < 2; ++y)  // y<4
                for (int x = 0; x < lx; ++x)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            if ((y < 1 && TStencil::Ci_y[idx] > 0) || (y >= 1 && TStencil::Ci_y[idx] < 0)) {
                                int k = 2 * lx * lz + z + x * lz + y * lx * lz;  // 2*lx*lz +
                                int yOffset = (y < 1) ? lw : ly - lw - 2;        // (y<2) ? lw-1 : ly-lw-3
                                int kGlobal = z + (y + yOffset) * lz + x * ly * lz;
                                obj.getDistribution()[kGlobal * TStencil::Q + idx] =
                                    obj.getCommDistribution()[k * TStencil::Q + idx];
                            }
                        }
        }
        if (TLattice::HaloZWidth) {
            lw = TLattice::HaloZWidth;
            for (int z = 0; z < 2; ++z)  // z<4
                for (int x = 0; x < lx; ++x)
                    for (int y = 0; y < ly; ++y)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            if ((z < 1 && TStencil::Ci_z[idx] > 0) || (z >= 1 && TStencil::Ci_z[idx] < 0)) {
                                int k = 2 * lx * ly + y + x * ly + z * lx * ly;  // 2*lx*ly +
                                int zOffset = (z < 1) ? lw : lz - lw - 2;        // (z<2) ? lw-1 : lz-lw-3
                                int kGlobal = (z + zOffset) + y * lz + x * ly * lz;
                                obj.getDistribution()[kGlobal * TStencil::Q + idx] =
                                    obj.getCommDistribution()[k * TStencil::Q + idx];
                            }
                        }
        }
    }
}

template <class TDerived, int TWidth>
template <class TLattice, class TDistribution>
inline void Parallel<TDerived, TWidth>::updateDistributionBeforeCommunicationAll(TDistribution& obj) {
#pragma omp master
    {
        using TStencil = typename TDistribution::Stencil;
        int lz = TLattice::LZdiv;
        int ly = TLattice::LYdiv;
        int lx = TLattice::LXdiv;
        int lw;
        if (TLattice::HaloXWidth) {
            lw = TLattice::HaloXWidth;
            for (int x = 0; x < 2 * lw; ++x)  // 4*lw
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = z + y * lz + x * ly * lz;
                            int xOffset = (x < lw) ? lw : lx - 3 * lw;  // (x<2*lw) ? 0 : lx-4*lw
                            int kGlobal = z + y * lz + (x + xOffset) * ly * lz;
                            // if (x<lw&&mpi.rank==4&&TIME==1000) std::cerr<<"before "<<kGlobal<<"
                            // "<<obj.getDistribution()[kGlobal * TStencil::Q + idx]<<std::endl;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getDistribution()[kGlobal * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloYWidth) {
            lw = TLattice::HaloYWidth;
            for (int y = 0; y < 2 * lw; ++y)  // 4*lw
                for (int x = 0; x < lx; ++x)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = z + x * lz + y * lx * lz;
                            int yOffset = (y < lw) ? lw : ly - 3 * lw;  // (y<2*lw) ? 0 : ly-4*lw
                            int kGlobal = z + (y + yOffset) * lz + x * ly * lz;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getDistribution()[kGlobal * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloZWidth) {
            lw = TLattice::HaloZWidth;
            for (int z = 0; z < 2 * lw; ++z)  // 4*lw
                for (int x = 0; x < lx; ++x)
                    for (int y = 0; y < ly; ++y)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = y + x * ly + z * lx * ly;
                            int zOffset = (z < lw) ? lw : lz - 3 * lw;  // (z<2*lw) ? 0 : lz-4*lw
                            int kGlobal = (z + zOffset) + y * lz + x * ly * lz;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getDistribution()[kGlobal * TStencil::Q + idx];
                        }
        }
    }
    return;
}

template <class TDerived, int TWidth>
template <class TLattice, class TDistribution>
inline void Parallel<TDerived, TWidth>::updateDistributionAfterCommunicationAll(TDistribution& obj) {
#pragma omp master
    {
        using TStencil = typename TDistribution::Stencil;
        int lz = TLattice::LZdiv;
        int ly = TLattice::LYdiv;
        int lx = TLattice::LXdiv;
        int lw;
        if (TLattice::HaloXWidth) {
            lw = TLattice::HaloXWidth;
            for (int x = 0; x < 2 * lw; ++x)  // 4*lw
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = 2 * lw * ly * lz + z + y * lz + x * ly * lz;  // 2*lw*ly*lz +
                            int xOffset = (x < lw) ? 0 : lx - 2 * lw;             // x<2*lw
                            int kGlobal = z + y * lz + (x + xOffset) * ly * lz;
                            // if (x>=lw&&mpi.rank==3&&TIME==1000) std::cerr<<"after "<<kGlobal<<"
                            // "<<obj.getCommDistribution()[k * TStencil::Q + idx]<<std::endl;
                            obj.getDistribution()[kGlobal * TStencil::Q + idx] =
                                obj.getCommDistribution()[k * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloYWidth) {
            lw = TLattice::HaloYWidth;
            for (int y = 0; y < 4 * lw; ++y)
                for (int x = 0; x < lx; ++x)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = 2 * lw * lx * lz + z + x * lz + y * lx * lz;  // 2*lw*lx*lz +
                            int yOffset = (y < lw) ? 0 : ly - 2 * lw;             // y<2*lw
                            int kGlobal = z + (y + yOffset) * lz + x * ly * lz;
                            obj.getDistribution()[kGlobal * TStencil::Q + idx] =
                                obj.getCommDistribution()[k * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloZWidth) {
            lw = TLattice::HaloZWidth;
            for (int z = 0; z < 4 * lw; ++z)
                for (int x = 0; x < lx; ++x)
                    for (int y = 0; y < ly; ++y)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = 2 * lw * lx * ly + y + x * ly + z * lx * ly;  // 2*lw*lx*ly +
                            int zOffset = (z < lw) ? 0 : lz - 2 * lw;             // z<2*lw
                            int kGlobal = (z + zOffset) + y * lz + x * ly * lz;
                            obj.getDistribution()[kGlobal * TStencil::Q + idx] =
                                obj.getCommDistribution()[k * TStencil::Q + idx];
                        }
        }
    }
    return;
}

template <class TDerived, int TWidth>
template <class TLattice, class TDistribution>
inline void Parallel<TDerived, TWidth>::updateDistributionBeforeCommunicationAllOld(TDistribution& obj) {
#pragma omp master
    {
        using TStencil = typename TDistribution::Stencil;
        int lz = TLattice::LZdiv;
        int ly = TLattice::LYdiv;
        int lx = TLattice::LXdiv;
        int lw;
        if (TLattice::HaloXWidth) {
            lw = TLattice::HaloXWidth;
            for (int x = 0; x < 2 * lw; ++x)  // 4*lw
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = z + y * lz + x * ly * lz;
                            int xOffset = (x < lw) ? lw : lx - 3 * lw;  // (x<2*lw) ? 0 : lx-4*lw
                            int kGlobal = z + y * lz + (x + xOffset) * ly * lz;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getDistributionOld()[kGlobal * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloYWidth) {
            lw = TLattice::HaloYWidth;
            for (int y = 0; y < 2 * lw; ++y)  // 4*lw
                for (int x = 0; x < lx; ++x)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = z + x * lz + y * lx * lz;
                            int yOffset = (y < lw) ? lw : ly - 3 * lw;  // (y<2*lw) ? 0 : ly-4*lw
                            int kGlobal = z + (y + yOffset) * lz + x * ly * lz;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getDistributionOld()[kGlobal * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloZWidth) {
            lw = TLattice::HaloZWidth;
            for (int z = 0; z < 2 * lw; ++z)  // 4*lw
                for (int x = 0; x < lx; ++x)
                    for (int y = 0; y < ly; ++y)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = y + x * ly + z * lx * ly;
                            int zOffset = (z < lw) ? lw : lz - 3 * lw;  // (z<2*lw) ? 0 : lz-4*lw
                            int kGlobal = (z + zOffset) + y * lz + x * ly * lz;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getDistributionOld()[kGlobal * TStencil::Q + idx];
                        }
        }
    }
    return;
}

template <class TDerived, int TWidth>
template <class TLattice, class TDistribution>
inline void Parallel<TDerived, TWidth>::updateDistributionAfterCommunicationAllOld(TDistribution& obj) {
#pragma omp master
    {
        using TStencil = typename TDistribution::Stencil;
        int lz = TLattice::LZdiv;
        int ly = TLattice::LYdiv;
        int lx = TLattice::LXdiv;
        int lw;
        if (TLattice::HaloXWidth) {
            lw = TLattice::HaloXWidth;
            for (int x = 0; x < 2 * lw; ++x)  // 4*lw
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = 2 * lw * ly * lz + z + y * lz + x * ly * lz;  // 2*lw*ly*lz +
                            int xOffset = (x < lw) ? 0 : lx - 4 * lw;             // x<2*lw
                            int kGlobal = z + y * lz + (x + xOffset) * ly * lz;
                            obj.getDistributionOld()[kGlobal * TStencil::Q + idx] =
                                obj.getCommDistribution()[k * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloYWidth) {
            lw = TLattice::HaloYWidth;
            for (int y = 0; y < 4 * lw; ++y)
                for (int x = 0; x < lx; ++x)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = 2 * lw * lx * lz + z + x * lz + y * lx * lz;  // 2*lw*lx*lz +
                            int yOffset = (y < lw) ? 0 : ly - 4 * lw;             // y<2*lw
                            int kGlobal = z + (y + yOffset) * lz + x * ly * lz;
                            obj.getDistributionOld()[kGlobal * TStencil::Q + idx] =
                                obj.getCommDistribution()[k * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloZWidth) {
            lw = TLattice::HaloZWidth;
            for (int z = 0; z < 4 * lw; ++z)
                for (int x = 0; x < lx; ++x)
                    for (int y = 0; y < ly; ++y)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = 2 * lw * lx * ly + y + x * ly + z * lx * ly;  // 2*lw*lx*ly +
                            int zOffset = (z < lw) ? 0 : lz - 4 * lw;             // z<2*lw
                            int kGlobal = (z + zOffset) + y * lz + x * ly * lz;
                            obj.getDistributionOld()[kGlobal * TStencil::Q + idx] =
                                obj.getCommDistribution()[k * TStencil::Q + idx];
                        }
        }
    }
    return;
}

template <class TDerived, int TWidth>
template <class TLattice, class TDistribution>
inline void Parallel<TDerived, TWidth>::updateDistributionBeforeCommunicationAllEquilibrium(TDistribution& obj) {
#pragma omp master
    {
        using TStencil = typename TDistribution::Stencil;
        int lz = TLattice::LZdiv;
        int ly = TLattice::LYdiv;
        int lx = TLattice::LXdiv;
        int lw;
        if (TLattice::HaloXWidth) {
            lw = TLattice::HaloXWidth;
            for (int x = 0; x < 2 * lw; ++x)  // 4*lw
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = z + y * lz + x * ly * lz;
                            int xOffset = (x < lw) ? lw : lx - 3 * lw;  // (x<2*lw) ? 0 : lx-4*lw
                            int kGlobal = z + y * lz + (x + xOffset) * ly * lz;
                            // std::cerr<<obj.getCommDistribution().size()<<" "<<k*TStencil::Q + idx<<std::endl;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getEquilibrium()[kGlobal * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloYWidth) {
            lw = TLattice::HaloYWidth;
            for (int y = 0; y < 2 * lw; ++y)  // 4*lw
                for (int x = 0; x < lx; ++x)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = z + x * lz + y * lx * lz;
                            int yOffset = (y < lw) ? lw : ly - 3 * lw;  // (y<2*lw) ? 0 : ly-4*lw
                            int kGlobal = z + (y + yOffset) * lz + x * ly * lz;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getEquilibrium()[kGlobal * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloZWidth) {
            lw = TLattice::HaloZWidth;
            for (int z = 0; z < 2 * lw; ++z)
                for (int x = 0; x < lx; ++x)
                    for (int y = 0; y < ly; ++y)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = y + x * ly + z * lx * ly;
                            int zOffset = (z < lw) ? lw : lz - 3 * lw;  // // (z<2*lw) ? 0 : lz-4*lw
                            int kGlobal = (z + zOffset) + y * lz + x * ly * lz;
                            obj.getCommDistribution()[k * TStencil::Q + idx] =
                                obj.getEquilibrium()[kGlobal * TStencil::Q + idx];
                        }
        }
    }
    return;
}

template <class TDerived, int TWidth>
template <class TLattice, class TDistribution>
inline void Parallel<TDerived, TWidth>::updateDistributionAfterCommunicationAllEquilibrium(TDistribution& obj) {
#pragma omp master
    {
        using TStencil = typename TDistribution::Stencil;
        int lz = TLattice::LZdiv;
        int ly = TLattice::LYdiv;
        int lx = TLattice::LXdiv;
        int lw;
        if (TLattice::HaloXWidth) {
            lw = TLattice::HaloXWidth;
            for (int x = 0; x < 2 * lw; ++x)  // 4*lw
                for (int y = 0; y < ly; ++y)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = 2 * lw * ly * lz + z + y * lz + x * ly * lz;  // 2*lw*ly*lz +
                            int xOffset = (x < lw) ? 0 : lx - 2 * lw;             // x<2*lw
                            int kGlobal = z + y * lz + (x + xOffset) * ly * lz;
                            obj.getEquilibrium()[kGlobal * TStencil::Q + idx] =
                                obj.getCommDistribution()[k * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloYWidth) {
            lw = TLattice::HaloYWidth;
            for (int y = 0; y < 4 * lw; ++y)
                for (int x = 0; x < lx; ++x)
                    for (int z = 0; z < lz; ++z)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = 2 * lw * lx * lz + z + x * lz + y * lx * lz;  // 2*lw*lx*lz +
                            int yOffset = (y < lw) ? 0 : ly - 2 * lw;             // y<2*lw
                            int kGlobal = z + (y + yOffset) * lz + x * ly * lz;
                            obj.getEquilibrium()[kGlobal * TStencil::Q + idx] =
                                obj.getCommDistribution()[k * TStencil::Q + idx];
                        }
        }
        if (TLattice::HaloZWidth) {
            lw = TLattice::HaloZWidth;
            for (int z = 0; z < 4 * lw; ++z)
                for (int x = 0; x < lx; ++x)
                    for (int y = 0; y < ly; ++y)
                        for (int idx = 0; idx < TStencil::Q; ++idx) {
                            int k = 2 * lw * lx * ly + y + x * ly + z * lx * ly;  // 2*lw*lx*ly +
                            int zOffset = (z < lw) ? 0 : lz - 2 * lw;             // z<2*lw
                            int kGlobal = (z + zOffset) + y * lz + x * ly * lz;
                            obj.getEquilibrium()[kGlobal * TStencil::Q + idx] =
                                obj.getCommDistribution()[k * TStencil::Q + idx];
                        }
        }
    }
    return;
}

#include "NoParallel.hh"
#include "ParallelX.hh"
#include "ParallelY.hh"
#include "ParallelZ.hh"

// TODO: Move child classes such ParallelX, ParallelXY, ParallelXYZ etc. into separate files and include them here
// #include "ParallelXY.hh"
// #include "ParallelXZ.hh"
// #include "ParallelYZ.hh"
// #include "ParallelXYZ.hh"
// TODO: Extend tests from ParallelX to other parallel patterns to ensure correctness of parallelisation
// The ultimate goal is to implement full 2D and 3D parallelisation
