#pragma once
#include <stdexcept>
#include "Mpi.hh"
#include "Global.hh"
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
template<class TDerived, int TNumNeighbors=1>
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
        template<class TLattice, class TParameter>
        inline void communicateParameter(TParameter& obj);

        /**
         * \brief Function to update the vector containing the parameters in the communication region before communication (implemented along X and Y).
         * \param obj Object of chosen parameter.
         */
        template<class TLattice, class TParameter>
        inline void updateParameterBeforeCommunication(TParameter& obj);

        /**
         * \brief Function to update the vector containing the parameters in the communication region after communication (implemented along X and Y).
         * \param obj Object of chosen parameter.
         */
        template<class TLattice, class TParameter>
        inline void updateParameterAfterCommunication(TParameter& obj);

        /**
         * \brief Function to update the vector containing the distribution in the communication region (implemented along X).
         * \param obj Object of the distribution.
         */
        template<class TLattice, class TDistribution>
        inline void communicateDistribution(TDistribution& obj);

        /**
         * \brief Function to update unknown distributions in the adjacent processors streamed from the edge.
         * \param obj Object of the distribution.
         */
        template<class TLattice, class TDistribution>
        inline void updateDistributionBeforeCommunication(TDistribution& obj);

        /**
         * \brief Function to update the vector containing the distribution for the communication region (implemented along X).
         * \param obj Object of the distribution.
         */
        template<class TLattice, class TDistribution>
        inline void updateDistributionAfterCommunication(TDistribution& obj);

        template<class TLattice, class TDistribution>
        inline void updateDistributionBeforeCommunicationAll(TDistribution& obj);

        template<class TLattice, class TDistribution>
        inline void updateDistributionAfterCommunicationAll(TDistribution& obj);

        template<class TLattice, class TDistribution>
        inline void updateDistributionBeforeCommunicationAllOld(TDistribution& obj);

        template<class TLattice, class TDistribution>
        inline void updateDistributionAfterCommunicationAllOld(TDistribution& obj);

        template<class TLattice, class TDistribution>
        inline void updateDistributionBeforeCommunicationAllEquilibrium(TDistribution& obj);

        template<class TLattice, class TDistribution>
        inline void updateDistributionAfterCommunicationAllEquilibrium(TDistribution& obj);

        template<class TLattice, class TDistribution>
        inline void communicateDistributionAll(TDistribution& obj);

        //template<class TLattice, class TDistribution>
        //inline void communicateDistributionAllOld(TDistribution& obj);

        static constexpr int NumNeighbors = TNumNeighbors;

    protected:
        int mMaxNeighbors = 0;
        // int mNumDirections = 0; //!<Number of communication channels. Must be redefined in child classes.
        std::vector<int> mNeighbors; //!<IDs of the neighboring processes.
        std::vector<int> mI0Send; //!<Lattice index to begin sending to each processor.
        std::vector<int> mI0Recv; //!<Lattice index to begin receiving from each processor.
        std::vector<int> mI0SendDistr; //!<Lattice index to begin sending to each processor for the distributions.
        std::vector<int> mI0RecvDistr; //!<Lattice index to begin receiving from each processor for the distributions.
};


template<class TDerived, int TNumNeighbors>
Parallel<TDerived,TNumNeighbors>::Parallel() {
    if(mMaxNeighbors < TNumNeighbors) mMaxNeighbors = TNumNeighbors;
}


/**
 * \details This will communicate the chosen parameter using MPI_Isend and MPI_Irecv, which are non-blocking methods of
 *          communication. This means that each process does not need to wait for the other processes to communicate. At
 *          the end of this function, we have a MPI_Waitall call, to ensure all processes are synced.
 */
template<class TDerived, int TNumNeighbors>
template<class TLattice, class TParameter>
inline void Parallel<TDerived,TNumNeighbors>::communicateParameter(TParameter& obj) {
    #ifdef MPIPARALLEL
    if (mpi.size == 1) return;
    //std::cout<<TParameter::mName<<" "<<obj.getCommParameter()[0]<<std::endl;
    //exit(1);
    #pragma omp master
    {
    int nNeighbors = mNeighbors.size();
    MPI_Request commrequest[2*nNeighbors];
    MPI_Status commstatus[2*nNeighbors];

    for (int iNeighbor=0; iNeighbor<nNeighbors; iNeighbor++) {

        int tag = iNeighbor;
        MPI_Isend(&obj.getCommParameter()[mI0Send[iNeighbor]*obj.mNum],
                  TNumNeighbors * TLattice::Face[iNeighbor/2] * obj.mNum,
                  mpi_get_type<typename TParameter::ParamType,TLattice>(),
                  mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2*iNeighbor]);

        tag = (iNeighbor%2==0) ? iNeighbor+1 : iNeighbor-1;
        MPI_Irecv(&obj.getCommParameter()[mI0Recv[iNeighbor]*obj.mNum],
                  TNumNeighbors * TLattice::Face[iNeighbor/2] * obj.mNum,
                  mpi_get_type<typename TParameter::ParamType,TLattice>(),
                  mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2*iNeighbor+1]);
        
    }

    MPI_Waitall(2*nNeighbors, commrequest, commstatus);
    }
    #endif
}

/*
template<class TDerived, int TNumNeighbors>
template<class TLattice, class TParameter>
inline void Parallel<TDerived,TNumNeighbors>::updateParameterBeforeCommunication(TParameter& obj) {
#pragma omp master
{
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int ln = 0;                   // number of elements in the communication region already handled
    if (TLattice::HaloXWidth) {
        int lw = TLattice::HaloXWidth;
        for(int x=0; x<2*lw; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                {
                    int k = z + y*lz + x*ly*lz;
                    int kComm = ln+k;
                    int xOffset = (x<lw) ? lw : lx-3*lw;
                    int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 2*lw*ly*lz;
    }
    if (TLattice::HaloYWidth) {
        int lw = TLattice::HaloYWidth;
        for(int x=0; x<lx; ++x)
            for(int y=0; y<2*lw; ++y)
                for(int z=0; z<lz; ++z)
                {
                    int k = z + x*lz + y*lx*lz;
                    int kComm = ln+k;
                    int yOffset = (y<lw) ? lw : ly-3*lw;
                    int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 2*lw*lx*lz;
    }
    if (TLattice::HaloZWidth) {
        int lw = TLattice::HaloZWidth;
        for(int x=0; x<lx; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<2*lw; ++z)
                {
                    int k = y + x*ly + z*lx*ly;
                    int kComm = ln+k;
                    int zOffset = (z<lw) ? lw : lz-3*lw;
                    int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 2*lw*lx*ly;
    }
}
}

template<class TDerived, int TNumNeighbors>
template<class TLattice, class TParameter>
inline void Parallel<TDerived,TNumNeighbors>::updateParameterAfterCommunication(TParameter& obj) {
#pragma omp master
{
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int ln = 0;                   // number of elements in the communication region already handled
    if (TLattice::HaloXWidth) {
        int lw = TLattice::HaloXWidth;
        for(int x=0; x<2*lw; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                {
                    int k = z + y*lz + x*ly*lz;
                    int kComm = 2*lw*ly*lz+ln+k;
                    int xOffset = (x<lw) ? 0 : lx-2*lw;
                    int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 2*lw*ly*lz;
    }
    if (TLattice::HaloYWidth) {
        int lw = TLattice::HaloYWidth;
        for(int x=0; x<lx; ++x)
            for(int y=0; y<2*lw; ++y)
                for(int z=0; z<lz; ++z)
                {
                    int k = z + x*lz + y*lx*lz;
                    int kComm = 2*lw*lx*lz+ln+k;
                    int yOffset = (y<lw) ? 0 : ly-2*lw;
                    int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 2*lw*lx*lz;
    }
    if (TLattice::HaloZWidth) {
        int lw = TLattice::HaloZWidth;
        for(int x=0; x<lx; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<4*lw; ++z)
                {
                    int k = y + x*ly + z*lx*ly;
                    int kComm = 2*lw*lx*ly+ln+k;
                    int zOffset = (z<lw) ? 0 : lz-2*lw;
                    int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 2*lw*lx*ly;
    }
}
}
*/

template<class TDerived, int TWidth>
template<class TLattice, class TParameter>
inline void Parallel<TDerived,TWidth>::updateParameterBeforeCommunication(TParameter& obj) {
#pragma omp master
{
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int ln = 0;                   // number of elements in the communication region already handled
    if (TLattice::HaloXWidth) {
        int lw = TLattice::HaloXWidth;
        for(int x=0; x<4*lw; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                {
                    int k = z + y*lz + x*ly*lz;
                    int kComm = ln+k;
                    int xOffset = (x<2*lw) ? 0 : lx-4*lw;
                    int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 4*lw*ly*lz;
    }
    if (TLattice::HaloYWidth) {
        int lw = TLattice::HaloYWidth;
        for(int x=0; x<lx; ++x)
            for(int y=0; y<4*lw; ++y)
                for(int z=0; z<lz; ++z)
                {
                    int k = z + x*lz + y*lx*lz;
                    int kComm = ln+k;
                    int yOffset = (y<2*lw) ? 0 : ly-4*lw;
                    int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 4*lw*lx*lz;
    }
    if (TLattice::HaloZWidth) {
        int lw = TLattice::HaloZWidth;
        for(int x=0; x<lx; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<4*lw; ++z)
                {
                    int k = y + x*ly + z*lx*ly;
                    int kComm = ln+k;
                    int zOffset = (z<2*lw) ? 0 : lz-4*lw;
                    int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 4*lw*lx*ly;
    }
}
}

template<class TDerived, int TWidth>
template<class TLattice, class TParameter>
inline void Parallel<TDerived,TWidth>::updateParameterAfterCommunication(TParameter& obj) {
#pragma omp master
{
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int ln = 0;                   // number of elements in the communication region already handled
    if (TLattice::HaloXWidth) {
        int lw = TLattice::HaloXWidth;
        for(int x=0; x<4*lw; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                {
                    int k = z + y*lz + x*ly*lz;
                    int kComm = ln+k;
                    int xOffset = (x<2*lw) ? 0 : lx-4*lw;
                    int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 4*lw*ly*lz;
    }
    if (TLattice::HaloYWidth) {
        int lw = TLattice::HaloYWidth;
        for(int x=0; x<lx; ++x)
            for(int y=0; y<4*lw; ++y)
                for(int z=0; z<lz; ++z)
                {
                    int k = z + x*lz + y*lx*lz;
                    int kComm = ln+k;
                    int yOffset = (y<2*lw) ? 0 : ly-4*lw;
                    int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 4*lw*lx*lz;
    }
    if (TLattice::HaloZWidth) {
        int lw = TLattice::HaloZWidth;
        for(int x=0; x<lx; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<4*lw; ++z)
                {
                    int k = y + x*ly + z*lx*ly;
                    int kComm = ln+k;
                    int zOffset = (z<2*lw) ? 0 : lz-4*lw;
                    int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                    for(int instance = 0; instance < TParameter::mInstances; instance++)
                        for(int direction = 0; direction < TParameter::mDirections; direction++)
                            obj.getParameter()[kGlobal*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction] =
                                obj.getCommParameter()[kComm*TParameter::mInstances*TParameter::mNum + (TParameter::mInstances>1)*instance*TParameter::mNum + (TParameter::mNum>1)*direction];
                }
        ln += 4*lw*lx*ly;
    }
}
}

template<class TDerived, int TNumNeighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,TNumNeighbors>::communicateDistribution(TDistribution& obj) {

    #ifdef MPIPARALLEL
    if (mpi.size == 1) return;

    using TStencil = typename TDistribution::Stencil;
    static MPI_Datatype distributionType = TDerived::template createDistributionType<TLattice,TStencil>();

    #pragma omp master
    {
    MPI_Request commdist_request[20];

    int id = 0;

    for(int dir = 0; dir < TDerived::mNumDirections; dir++)
    
        for(int idx = 1; idx < TStencil::Q; idx++) {

            int direction = TStencil::Ci_xyz(TDerived::mCommDirection[dir])[idx]; // 0 --> TDerived::mCommDirection[dir]
            int iNeighbor = 2*dir;
            int tag = idx; // iNeighbor;

            if(direction == -1) {
                MPI_Isend(&obj.getCommDistribution()[mI0SendDistr[0]*TStencil::Q + idx],
                        1, distributionType, mNeighbors[iNeighbor],
                        tag, MPI_COMM_WORLD, &commdist_request[id++]);
                MPI_Irecv(&obj.getCommDistribution()[mI0RecvDistr[1]*TStencil::Q + idx],
                        1, distributionType, mNeighbors[iNeighbor+1],
                        tag, MPI_COMM_WORLD, &commdist_request[id++]);

            }
            else if(direction == 1) {
            
                MPI_Isend(&obj.getCommDistribution()[mI0SendDistr[1]*TStencil::Q + idx],
                        1, distributionType, mNeighbors[iNeighbor+1],
                        tag, MPI_COMM_WORLD, &commdist_request[id++]);
                MPI_Irecv(&obj.getCommDistribution()[mI0RecvDistr[0]*TStencil::Q + idx],
                        1, distributionType, mNeighbors[iNeighbor],
                        tag, MPI_COMM_WORLD, &commdist_request[id++]);

            }

        }

    MPI_Waitall(id, commdist_request, MPI_STATUSES_IGNORE);
    }
    #endif
}

template<class TDerived, int TNumNeighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,TNumNeighbors>::communicateDistributionAll(TDistribution& obj) {

    #ifdef MPIPARALLEL
    if (mpi.size == 1) return;

    using TStencil = typename TDistribution::Stencil;

    #pragma omp master
    {
    int nNeighbors = mNeighbors.size();
    MPI_Request commrequest[2*nNeighbors];

    for (int iNeighbor=0; iNeighbor<nNeighbors; iNeighbor++) {
        int tag = iNeighbor;
        MPI_Isend(&obj.getCommDistribution()[mI0Send[iNeighbor]*TStencil::Q],
                  TNumNeighbors * TLattice::LY * TLattice::LZ * TStencil::Q,
                  mpi_get_type<double,TLattice>(),
                  mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2*iNeighbor]);

        tag = (iNeighbor%2==0) ? iNeighbor+1 : iNeighbor-1;
        MPI_Irecv(&obj.getCommDistribution()[mI0Recv[iNeighbor]*TStencil::Q],
                  TNumNeighbors * TLattice::LY * TLattice::LZ * TStencil::Q,
                  mpi_get_type<double,TLattice>(),
                  mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2*iNeighbor+1]);
    }

    MPI_Waitall(2*nNeighbors, commrequest, MPI_STATUSES_IGNORE);
    }
    #endif

}

template<class TDerived, int TNumNeighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,TNumNeighbors>::updateDistributionBeforeCommunication(TDistribution& obj) {
#pragma omp master
{
    using TStencil = typename TDistribution::Stencil;
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int lw;
    if (TLattice::HaloXWidth) {
        lw = TLattice::HaloXWidth;
        for(int x=0; x<2; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = z + y*lz + x*ly*lz;
                        int xOffset = (x<1) ? lw-1 : lx-lw-1;
                        int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                       // std::cerr<<obj.getCommDistribution().size()<<" "<<k*TStencil::Q + idx<<std::endl;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getDistribution()[kGlobal*TStencil::Q + idx];
                        //if (mpi.rank==2) std::cout<<k<<" "<<kGlobal<<" "<<obj.getDistribution()[kGlobal*TStencil::Q + idx]<<std::endl;
                        //if (mpi.rank==2) std::cout<<k<<" "<<obj.getCommDistribution()[k*TStencil::Q + idx]<<" "<<kGlobal<<" "<<obj.getDistribution()[kGlobal*TStencil::Q + idx]<<std::endl;
                    }
    }
    if (TLattice::HaloYWidth) {
        lw = TLattice::HaloYWidth;
        for(int y=0; y<2; ++y)
            for(int x=0; x<lx; ++x)
                for(int z=0; z<lz; ++z)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = z + x*lz + y*lx*lz;
                        int yOffset = (y<1) ? lw-1 : ly-lw-1;
                        int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getDistribution()[kGlobal*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloZWidth) {
        lw = TLattice::HaloZWidth;
        for(int z=0; z<2; ++z)
            for(int x=0; x<lx; ++x)
                for(int y=0; y<ly; ++y)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = y + x*ly + z*lx*ly;
                        int zOffset = (z<1) ? lw-1 : lz-lw-1;
                        int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getDistribution()[kGlobal*TStencil::Q + idx];
                    }
    }
}
}

template<class TDerived, int TNumNeighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,TNumNeighbors>::updateDistributionAfterCommunication(TDistribution& obj) {
#pragma omp master
{
    using TStencil = typename TDistribution::Stencil;
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int lw;
    if (TLattice::HaloXWidth) {
        lw = TLattice::HaloXWidth;
        for(int x=0; x<2; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                    for (int idx = 0; idx < TStencil::Q; ++idx)
                    {
                        if((x<1 && TStencil::Ci_x[idx]>0)||(x>=1 && TStencil::Ci_x[idx]<0)){
                            int k = 2*ly*lz + z + y*lz + x*ly*lz;
                            int xOffset = (x<1) ? lw : lx-lw-2;
                            int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                            obj.getDistribution()[kGlobal*TStencil::Q + idx] =
                                obj.getCommDistribution()[k*TStencil::Q + idx];
                            //if (mpi.rank==2) std::cout<<k<<" "<<obj.getCommDistribution()[k*TStencil::Q + idx]<<" "<<kGlobal<<" "<<obj.getDistribution()[kGlobal*TStencil::Q + idx]<<std::endl;
                        }
                        
                    }
    }
    if (TLattice::HaloYWidth) {
        lw = TLattice::HaloYWidth;
        for(int y=0; y<2; ++y)
            for(int x=0; x<lx; ++x)
                for(int z=0; z<lz; ++z)
                    for (int idx = 0; idx < TStencil::Q; ++idx)
                    {
                        if((y<1 && TStencil::Ci_y[idx]>0)||(y>=1 && TStencil::Ci_y[idx]<0)){
                            int k = 2*lx*lz + z + x*lz + y*lx*lz;
                            int yOffset = (y<1) ? lw : ly-lw-2;
                            int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                            obj.getDistribution()[kGlobal*TStencil::Q + idx] =
                                obj.getCommDistribution()[k*TStencil::Q + idx];
                        }
                    }
    }
    if (TLattice::HaloZWidth) {
        lw = TLattice::HaloZWidth;
        for(int z=0; z<2; ++z)
            for(int x=0; x<lx; ++x)
                for(int y=0; y<ly; ++y)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        if((z<1 && TStencil::Ci_z[idx]>0)||(z>=1 && TStencil::Ci_z[idx]<0)){
                            int k = 2*lx*ly + y + x*ly + z*lx*ly;
                            int zOffset = (z<1) ? lw : lz-lw-2;
                            int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                            obj.getDistribution()[kGlobal*TStencil::Q + idx] =
                                obj.getCommDistribution()[k*TStencil::Q + idx];
                        }
                    }
    }
}
}

template<class TDerived, int TNumNeighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,TNumNeighbors>::updateDistributionBeforeCommunicationAll(TDistribution& obj) {

    #pragma omp master
    {
    using TStencil = typename TDistribution::Stencil;
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int lw;
    if (TLattice::HaloXWidth) {
        lw = TLattice::HaloXWidth;
        for(int x=0; x<2*lw; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = z + y*lz + x*ly*lz;
                        int xOffset = (x<lw) ? lw : lx-3*lw;
                        int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                       // std::cerr<<obj.getCommDistribution().size()<<" "<<k*TStencil::Q + idx<<std::endl;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getDistribution()[kGlobal*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloYWidth) {
        lw = TLattice::HaloYWidth;
        for(int y=0; y<2*lw; ++y)
            for(int x=0; x<lx; ++x)
                for(int z=0; z<lz; ++z)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = z + x*lz + y*lx*lz;
                        int yOffset = (y<lw) ? lw : ly-3*lw;
                        int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getDistribution()[kGlobal*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloZWidth) {
        lw = TLattice::HaloZWidth;
        for(int z=0; z<2*lw; ++z)
            for(int x=0; x<lx; ++x)
                for(int y=0; y<ly; ++y)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = y + x*ly + z*lx*ly;
                        int zOffset = (z<lw) ? lw : lz-3*lw;
                        int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getDistribution()[kGlobal*TStencil::Q + idx];
                    }
    }
    }
    return;
    
}

template<class TDerived, int TNumNeighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,TNumNeighbors>::updateDistributionAfterCommunicationAll(TDistribution& obj) {
    #pragma omp master
    {
    using TStencil = typename TDistribution::Stencil;
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int lw;
    if (TLattice::HaloXWidth) {
        lw = TLattice::HaloXWidth;
        for(int x=0; x<2*lw; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                    for (int idx = 0; idx < TStencil::Q; ++idx)
                    {
                        int k = 2*lw*ly*lz + z + y*lz + x*ly*lz;
                        int xOffset = (x<lw) ? 0 : lx-4*lw;
                        int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                        obj.getDistribution()[kGlobal*TStencil::Q + idx] =
                            obj.getCommDistribution()[k*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloYWidth) {
        lw = TLattice::HaloYWidth;
        for(int y=0; y<4*lw; ++y)
            for(int x=0; x<lx; ++x)
                for(int z=0; z<lz; ++z)
                    for (int idx = 0; idx < TStencil::Q; ++idx)
                    {
                        int k = 2*lw*lx*lz + z + x*lz + y*lx*lz;
                        int yOffset = (y<lw) ? 0 : ly-4*lw;
                        int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                        obj.getDistribution()[kGlobal*TStencil::Q + idx] =
                            obj.getCommDistribution()[k*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloZWidth) {
        lw = TLattice::HaloZWidth;
        for(int z=0; z<4*lw; ++z)
            for(int x=0; x<lx; ++x)
                for(int y=0; y<ly; ++y)
                    for (int idx = 0; idx < TStencil::Q; ++idx)
                    {
                        int k = 2*lw*lx*ly + y + x*ly + z*lx*ly;
                        int zOffset = (z<lw) ? 0 : lz-4*lw;
                        int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                        obj.getDistribution()[kGlobal*TStencil::Q + idx] =
                            obj.getCommDistribution()[k*TStencil::Q + idx];
                    }
    }
    }
    return;
}


template<class TDerived, int TNumNeighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,TNumNeighbors>::updateDistributionBeforeCommunicationAllOld(TDistribution& obj) {

    #pragma omp master
    {
    using TStencil = typename TDistribution::Stencil;
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int lw;
    if (TLattice::HaloXWidth) {
        lw = TLattice::HaloXWidth;
        for(int x=0; x<2*lw; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = z + y*lz + x*ly*lz;
                        int xOffset = (x<lw) ? lw : lx-3*lw;
                        int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                       // std::cerr<<obj.getCommDistribution().size()<<" "<<k*TStencil::Q + idx<<std::endl;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getDistributionOld()[kGlobal*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloYWidth) {
        lw = TLattice::HaloYWidth;
        for(int y=0; y<2*lw; ++y)
            for(int x=0; x<lx; ++x)
                for(int z=0; z<lz; ++z)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = z + x*lz + y*lx*lz;
                        int yOffset = (y<lw) ? lw : ly-3*lw;
                        int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getDistributionOld()[kGlobal*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloZWidth) {
        lw = TLattice::HaloZWidth;
        for(int z=0; z<2*lw; ++z)
            for(int x=0; x<lx; ++x)
                for(int y=0; y<ly; ++y)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = y + x*ly + z*lx*ly;
                        int zOffset = (z<lw) ? lw : lz-3*lw;
                        int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getDistributionOld()[kGlobal*TStencil::Q + idx];
                    }
    }
    }
    return;
    
}

template<class TDerived, int TNumNeighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,TNumNeighbors>::updateDistributionAfterCommunicationAllOld(TDistribution& obj) {
    #pragma omp master
    {
    using TStencil = typename TDistribution::Stencil;
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int lw;
    if (TLattice::HaloXWidth) {
        lw = TLattice::HaloXWidth;
        for(int x=0; x<2*lw; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                    for (int idx = 0; idx < TStencil::Q; ++idx)
                    {
                        int k = 2*lw*ly*lz + z + y*lz + x*ly*lz;
                        int xOffset = (x<lw) ? 0 : lx-4*lw;
                        int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                        obj.getDistributionOld()[kGlobal*TStencil::Q + idx] =
                            obj.getCommDistribution()[k*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloYWidth) {
        lw = TLattice::HaloYWidth;
        for(int y=0; y<4*lw; ++y)
            for(int x=0; x<lx; ++x)
                for(int z=0; z<lz; ++z)
                    for (int idx = 0; idx < TStencil::Q; ++idx)
                    {
                        int k = 2*lw*lx*lz + z + x*lz + y*lx*lz;
                        int yOffset = (y<lw) ? 0 : ly-4*lw;
                        int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                        obj.getDistributionOld()[kGlobal*TStencil::Q + idx] =
                            obj.getCommDistribution()[k*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloZWidth) {
        lw = TLattice::HaloZWidth;
        for(int z=0; z<4*lw; ++z)
            for(int x=0; x<lx; ++x)
                for(int y=0; y<ly; ++y)
                    for (int idx = 0; idx < TStencil::Q; ++idx)
                    {
                        int k = 2*lw*lx*ly + y + x*ly + z*lx*ly;
                        int zOffset = (z<lw) ? 0 : lz-4*lw;
                        int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                        obj.getDistributionOld()[kGlobal*TStencil::Q + idx] =
                            obj.getCommDistribution()[k*TStencil::Q + idx];
                    }
    }
    }
    return;
}


template<class TDerived, int TNumNeighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,TNumNeighbors>::updateDistributionBeforeCommunicationAllEquilibrium(TDistribution& obj) {

    #pragma omp master
    {
    using TStencil = typename TDistribution::Stencil;
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int lw;
    if (TLattice::HaloXWidth) {
        lw = TLattice::HaloXWidth;
        for(int x=0; x<2*lw; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = z + y*lz + x*ly*lz;
                        int xOffset = (x<lw) ? lw : lx-3*lw;
                        int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                       // std::cerr<<obj.getCommDistribution().size()<<" "<<k*TStencil::Q + idx<<std::endl;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getEquilibrium()[kGlobal*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloYWidth) {
        lw = TLattice::HaloYWidth;
        for(int y=0; y<2*lw; ++y)
            for(int x=0; x<lx; ++x)
                for(int z=0; z<lz; ++z)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = z + x*lz + y*lx*lz;
                        int yOffset = (y<lw) ? lw : ly-3*lw;
                        int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getEquilibrium()[kGlobal*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloZWidth) {
        lw = TLattice::HaloZWidth;
        for(int z=0; z<2*lw; ++z)
            for(int x=0; x<lx; ++x)
                for(int y=0; y<ly; ++y)
                    for (int idx=0; idx<TStencil::Q; ++idx)
                    {
                        int k = y + x*ly + z*lx*ly;
                        int zOffset = (z<lw) ? lw : lz-3*lw;
                        int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                        obj.getCommDistribution()[k*TStencil::Q + idx] =
                            obj.getEquilibrium()[kGlobal*TStencil::Q + idx];
                    }
    }
    }
    return;
    
}

template<class TDerived, int TNumNeighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,TNumNeighbors>::updateDistributionAfterCommunicationAllEquilibrium(TDistribution& obj) {
    #pragma omp master
    {
    using TStencil = typename TDistribution::Stencil;
    int lz = TLattice::LZdiv;
    int ly = TLattice::LYdiv;
    int lx = TLattice::LXdiv;
    int lw;
    if (TLattice::HaloXWidth) {
        lw = TLattice::HaloXWidth;
        for(int x=0; x<2*lw; ++x)
            for(int y=0; y<ly; ++y)
                for(int z=0; z<lz; ++z)
                    for (int idx = 0; idx < TStencil::Q; ++idx)
                    {
                        int k = 2*lw*ly*lz + z + y*lz + x*ly*lz;
                        int xOffset = (x<lw) ? 0 : lx-4*lw;
                        int kGlobal = z + y*lz + (x+xOffset)*ly*lz;
                        obj.getEquilibrium()[kGlobal*TStencil::Q + idx] =
                            obj.getCommDistribution()[k*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloYWidth) {
        lw = TLattice::HaloYWidth;
        for(int y=0; y<4*lw; ++y)
            for(int x=0; x<lx; ++x)
                for(int z=0; z<lz; ++z)
                    for (int idx = 0; idx < TStencil::Q; ++idx)
                    {
                        int k = 2*lw*lx*lz + z + x*lz + y*lx*lz;
                        int yOffset = (y<lw) ? 0 : ly-4*lw;
                        int kGlobal = z + (y+yOffset)*lz + x*ly*lz;
                        obj.getEquilibrium()[kGlobal*TStencil::Q + idx] =
                            obj.getCommDistribution()[k*TStencil::Q + idx];
                    }
    }
    if (TLattice::HaloZWidth) {
        lw = TLattice::HaloZWidth;
        for(int z=0; z<4*lw; ++z)
            for(int x=0; x<lx; ++x)
                for(int y=0; y<ly; ++y)
                    for (int idx = 0; idx < TStencil::Q; ++idx)
                    {
                        int k = 2*lw*lx*ly + y + x*ly + z*lx*ly;
                        int zOffset = (z<lw) ? 0 : lz-4*lw;
                        int kGlobal = (z+zOffset) + y*lz + x*ly*lz;
                        obj.getEquilibrium()[kGlobal*TStencil::Q + idx] =
                            obj.getCommDistribution()[k*TStencil::Q + idx];
                    }
    }
    }
    return;
}


/**TDistribution&
 * \brief NoParallel is a dummy implementation of Parallel that does not split the lattice or perform any communication.
 */
class NoParallel : public Parallel<NoParallel,0> {
    public:
        // temporary implementation before a better idea
        static constexpr int mNumDirections = 0; // number of communication directions in 1D parallelisation
        static constexpr int mCommDirection[1] = {0}; // communication direction along X

        template<class TLattice>
        void init() {};

        #ifdef MPIPARALLEL
        template<class TLattice, class TStencil>
        static MPI_Datatype createDistributionType() { return MPI_DOUBLE; };
        #endif
};


/**
 * \brief ParallelX contains functions and data for MPI parallelisation divided evenly in the X direction.
 * This class contains functions to communicate parameters for gradient calculations and to communicate
 * distributions for streaming.
 */
template<int TNumNeighbors=1>
class ParallelX : public Parallel<ParallelX<TNumNeighbors>,TNumNeighbors> {
    public:
        
        static constexpr int mNumDirections = 1; // number of communication directions in 1D parallelisation
        static constexpr int mCommDirection[mNumDirections] = {0}; // communication direction along X

        /**
         * \brief Initialise MPI variables for the parallelisation method.
         */
        template<class TLattice>
        void init();

        #ifdef MPIPARALLEL
        template<class TLattice, class TStencil>
        static MPI_Datatype createDistributionType(); //!<Create datatype for streaming distributions
        #endif
};


template<int TNumNeighbors>
template<class TLattice>
void ParallelX<TNumNeighbors>::init() {
    if (mpi.size <= 1) return;

    const int faceSize = TLattice::Face[0] = TLattice::LY * TLattice::LZ;

    // Split lattice
    TLattice::HaloXWidth = TLattice::Neighbors = this->mMaxNeighbors;
    TLattice::HaloSize = faceSize * this->mMaxNeighbors;
    //std::cout<<this->mMaxNeighbors<<std::endl;
    if (TLattice::LX % mpi.size == 0) {
        TLattice::LXdiv = (TLattice::LX/mpi.size + 2*TNumNeighbors);
        TLattice::LXMPIOffset = (TLattice::LXdiv - 2*TNumNeighbors) * mpi.rank;
        TLattice::subArray[0] = TLattice::LX / mpi.size;
    }
    else{
        std::string err_message = "Currently, the size of the domain in the X direction must be divisible by the number of mpi ranks.";
        throw std::runtime_error(err_message);
    }
    TLattice::N = TLattice::LXdiv * faceSize;

    // Define neighbor processors
    int leftNeighbor = mpi.rank - 1;
    if (leftNeighbor == -1) leftNeighbor = mpi.size - 1;
    int rightNeighbor = mpi.rank + 1;
    if (rightNeighbor == mpi.size) rightNeighbor = 0;
    this->mNeighbors = {leftNeighbor, rightNeighbor};

    // Create communication objects
    this->mI0Send = std::vector<int>(2);
    this->mI0Send[0] = TNumNeighbors*faceSize;//TNumNeighbors*faceSize;
    this->mI0Send[1] = 2*TNumNeighbors*faceSize;//2*TNumNeighbors*faceSize;
    this->mI0Recv = std::vector<int>(2);
    this->mI0Recv[0] = 0;
    this->mI0Recv[1] = 3*TNumNeighbors*faceSize;

    this->mI0SendDistr = std::vector<int>(2);
    this->mI0SendDistr[0] = 0;
    this->mI0SendDistr[1] = faceSize;
    this->mI0RecvDistr = std::vector<int>(2);
    this->mI0RecvDistr[0] = 2*faceSize;
    this->mI0RecvDistr[1] = 3*faceSize;

    // Create MPI buffer
    #ifdef MPIPARALLEL
    const int bufSize = (faceSize * TNumNeighbors * (5 + 2) * 2 * 2 + 1000) * sizeof(double);

    if(bufSize > MPIBUFFERSIZE) {
        if(MPIBUFFERSIZE != 0) MPI_Buffer_detach(MPIBUFFER, &MPIBUFFERSIZE);
        if(MPIBUFFERSIZE != 0) delete[] MPIBUFFER;
        MPIBUFFER = new char[bufSize];
        MPI_Buffer_attach(MPIBUFFER, bufSize);
        MPIBUFFERSIZE = bufSize;
    }
    #endif
}


#ifdef MPIPARALLEL
template<int TNumNeighbors>
template<class TLattice, class TStencil>
MPI_Datatype ParallelX<TNumNeighbors>::createDistributionType() {
    const int faceSize = TLattice::LY*TLattice::LZ;
    MPI_Datatype distributionType;
    MPI_Type_vector(faceSize, 1, TStencil::Q, mpi_get_type<double,TLattice>(), &distributionType);
    MPI_Type_commit(&distributionType);
    return distributionType;
}
#endif


/**
 * \brief ParallelY contains functions and data for MPI parallelisation divided evenly in the Y direction.
 * This class contains functions to communicate parameters for gradient calculations and to communicate
 * distributions for streaming.
 */
template<int TNumNeighbors=1>
class ParallelY : public Parallel<ParallelY<TNumNeighbors>,TNumNeighbors> {
    public:

        static constexpr int mNumDirections = 1; // number of communication directions in 1D parallelisation
        static constexpr int mCommDirection[mNumDirections] = {1}; // communication direction along Y

        /**
         * \brief Initialise MPI variables for the parallelisation method.
         */
        template<class TLattice>
        void init();

        #ifdef MPIPARALLEL
        template<class TLattice, class TStencil>
        static MPI_Datatype createDistributionType(); //!<Create datatype for streaming distributions
        #endif
};


template<int TNumNeighbors>
template<class TLattice>
void ParallelY<TNumNeighbors>::init() {
    if (mpi.size <= 1) return;

    const int faceSize = TLattice::Face[0] = TLattice::LX * TLattice::LZ;

    // Split lattice
    TLattice::HaloYWidth = TLattice::Neighbors = this->mMaxNeighbors;
    TLattice::HaloSize = faceSize * this->mMaxNeighbors;

    if (TLattice::LY % mpi.size == 0) {
        TLattice::LYdiv = (TLattice::LY/mpi.size + 2*TNumNeighbors);
        TLattice::LYMPIOffset = (TLattice::LYdiv - 2*TNumNeighbors) * mpi.rank;
        TLattice::subArray[1] = TLattice::LY / mpi.size;
    }
    else{
        std::string err_message = "Currently, the size of the domain in the Y direction must be divisible by the number of mpi ranks.";
        throw std::runtime_error(err_message);
    }
    TLattice::N = TLattice::LYdiv * faceSize;

    // Define neighbor processors
    int bottomNeighbor = mpi.rank - 1;
    if (bottomNeighbor == -1) bottomNeighbor = mpi.size - 1;
    int topNeighbor = mpi.rank + 1;
    if (topNeighbor == mpi.size) topNeighbor = 0;
    this->mNeighbors = {bottomNeighbor, topNeighbor};

    // Create communication objects
    this->mI0Send = std::vector<int>(2);
    this->mI0Send[0] = TNumNeighbors*faceSize;//TNumNeighbors*faceSize;
    this->mI0Send[1] = 2*TNumNeighbors*faceSize;//2*TNumNeighbors*faceSize;
    this->mI0Recv = std::vector<int>(2);
    this->mI0Recv[0] = 0;
    this->mI0Recv[1] = 3*TNumNeighbors*faceSize;

    this->mI0SendDistr = std::vector<int>(2);
    this->mI0SendDistr[0] = 0;
    this->mI0SendDistr[1] = faceSize;
    this->mI0RecvDistr = std::vector<int>(2);
    this->mI0RecvDistr[0] = 2*faceSize;
    this->mI0RecvDistr[1] = 3*faceSize;

    // Create MPI buffer
    #ifdef MPIPARALLEL
    const int bufSize = (faceSize * TNumNeighbors * (5 + 2) * 2 * 2 + 1000) * sizeof(double);

    if(bufSize > MPIBUFFERSIZE) {
        if(MPIBUFFERSIZE != 0) MPI_Buffer_detach(MPIBUFFER, &MPIBUFFERSIZE);
        if(MPIBUFFERSIZE != 0) delete[] MPIBUFFER;
        MPIBUFFER = new char[bufSize];
        MPI_Buffer_attach(MPIBUFFER, bufSize);
        MPIBUFFERSIZE = bufSize;
    }
    #endif
}


#ifdef MPIPARALLEL
template<int TNumNeighbors>
template<class TLattice, class TStencil>
MPI_Datatype ParallelY<TNumNeighbors>::createDistributionType() {
    const int faceSize = TLattice::LX*TLattice::LZ;
    MPI_Datatype distributionType;
    MPI_Type_vector(faceSize, 1, TStencil::Q, mpi_get_type<double,TLattice>(), &distributionType);
    MPI_Type_commit(&distributionType);
    return distributionType;
}
#endif


/**
 * \brief ParallelZ contains functions and data for MPI parallelisation divided evenly in the Z direction.
 * This class contains functions to communicate parameters for gradient calculations and to communicate
 * distributions for streaming.
 */
template<int TNumNeighbors=1>
class ParallelZ : public Parallel<ParallelZ<TNumNeighbors>,TNumNeighbors> {
    public:

        static constexpr int mNumDirections = 1; // number of communication directions in 1D parallelisation
        static constexpr int mCommDirection[mNumDirections] = {2}; // communication direction along Z

        /**
         * \brief Initialise MPI variables for the parallelisation method.
         */
        template<class TLattice>
        void init();

        #ifdef MPIPARALLEL
        template<class TLattice, class TStencil>
        static MPI_Datatype createDistributionType(); //!<Create datatype for streaming distributions
        #endif
};


template<int TNumNeighbors>
template<class TLattice>
void ParallelZ<TNumNeighbors>::init() {
    if (mpi.size <= 1) return;

    const int faceSize = TLattice::Face[0] = TLattice::LX * TLattice::LY;

    // Split lattice
    TLattice::HaloZWidth = TLattice::Neighbors = this->mMaxNeighbors;
    TLattice::HaloSize = faceSize * this->mMaxNeighbors;

    if (TLattice::LZ % mpi.size == 0) {
        TLattice::LZdiv = (TLattice::LZ/mpi.size + 2*TNumNeighbors);
        TLattice::LZMPIOffset = (TLattice::LZdiv - 2*TNumNeighbors) * mpi.rank;
        TLattice::subArray[2] = TLattice::LZ / mpi.size;
    }
    else{
        std::string err_message = "Currently, the size of the domain in the Z direction must be divisible by the number of mpi ranks.";
        throw std::runtime_error(err_message);
    }
    TLattice::N = TLattice::LZdiv * faceSize;

    // Define neighbor processors
    int frontNeighbor = mpi.rank - 1;
    if (frontNeighbor == -1) frontNeighbor = mpi.size - 1;
    int backNeighbor = mpi.rank + 1;
    if (backNeighbor == mpi.size) backNeighbor = 0;
    this->mNeighbors = {frontNeighbor, backNeighbor};

    // Create communication objects
    this->mI0Send = std::vector<int>(2);
    this->mI0Send[0] = TNumNeighbors*faceSize;//TNumNeighbors*faceSize;
    this->mI0Send[1] = 2*TNumNeighbors*faceSize;//2*TNumNeighbors*faceSize;
    this->mI0Recv = std::vector<int>(2);
    this->mI0Recv[0] = 0;
    this->mI0Recv[1] = 3*TNumNeighbors*faceSize;

    this->mI0SendDistr = std::vector<int>(2);
    this->mI0SendDistr[0] = 0;
    this->mI0SendDistr[1] = faceSize;
    this->mI0RecvDistr = std::vector<int>(2);
    this->mI0RecvDistr[0] = 2*faceSize;
    this->mI0RecvDistr[1] = 3*faceSize;

    // Create MPI buffer
    #ifdef MPIPARALLEL
    const int bufSize = (faceSize * TNumNeighbors * (5 + 2) * 2 * 2 + 1000) * sizeof(double);

    if(bufSize > MPIBUFFERSIZE) {
        if(MPIBUFFERSIZE != 0) MPI_Buffer_detach(MPIBUFFER, &MPIBUFFERSIZE);
        if(MPIBUFFERSIZE != 0) delete[] MPIBUFFER;
        MPIBUFFER = new char[bufSize];
        MPI_Buffer_attach(MPIBUFFER, bufSize);
        MPIBUFFERSIZE = bufSize;
    }
    #endif
}


#ifdef MPIPARALLEL
template<int TNumNeighbors>
template<class TLattice, class TStencil>
MPI_Datatype ParallelZ<TNumNeighbors>::createDistributionType() {
    const int faceSize = TLattice::LX*TLattice::LY;
    MPI_Datatype distributionType;
    MPI_Type_vector(faceSize, 1, TStencil::Q, mpi_get_type<double,TLattice>(), &distributionType);
    MPI_Type_commit(&distributionType);
    return distributionType;
}
#endif


// TODO: Move child classes such ParallelX, ParallelXY, ParallelXYZ etc. into separate files and include them here
// #include "ParallelX.hh"
// #include "ParallelY.hh"
// #include "ParallelZ.hh"
// #include "ParallelXY.hh"
// #include "ParallelXZ.hh"
// #include "ParallelYZ.hh"
// #include "ParallelXYZ.hh"
// TODO: Extend tests from ParallelX to other parallel patterns to ensure correctness of parallelisation
// The ultimate goal is to implement full 2D and 3D parallelisation