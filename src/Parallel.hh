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
         *        from X_Parallel to here.
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
        inline void communicate(TParameter& obj);

        /**
         * \brief Function to update unknown distributions in the adjacent processors streamed from the edge.
         * \param obj Object of the distribution.
         */
        template<class TLattice, class TDistribution>
        inline void communicateDistribution(TDistribution& obj);

        template<class TLattice, class TDistribution>
        inline void communicateDistributionAll(TDistribution& obj);

        template<class TLattice, class TDistribution>
        inline void communicateDistributionAllOld(TDistribution& obj);

        static constexpr int NumNeighbors = TNumNeighbors;

    protected:
        int mMaxNeighbors = 0;
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
inline void Parallel<TDerived,TNumNeighbors>::communicate(TParameter& obj) {
    #ifdef MPIPARALLEL
    if (mpi.size == 1) return;

    #pragma omp master
    {
    int nNeighbors = mNeighbors.size();
    MPI_Request commrequest[2*nNeighbors];

    for (int iNeighbor=0; iNeighbor<nNeighbors; iNeighbor++) {
        int tag = iNeighbor;
        MPI_Isend(&obj.mv_Parameter[mI0Send[iNeighbor]*obj.mNum],
                  TNumNeighbors * TLattice::LY * TLattice::LZ * obj.mNum,
                  mpi_get_type<typename TParameter::ParamType>(),
                  mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2*iNeighbor]);

        tag = (iNeighbor%2==0) ? iNeighbor+1 : iNeighbor-1;
        MPI_Irecv(&obj.mv_Parameter[mI0Recv[iNeighbor]*obj.mNum],
                  TNumNeighbors * TLattice::LY * TLattice::LZ * obj.mNum,
                  mpi_get_type<typename TParameter::ParamType>(),
                  mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2*iNeighbor+1]);
    }

    MPI_Waitall(2*nNeighbors, commrequest, MPI_STATUSES_IGNORE);
    }
    #endif
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

    int leftorright;
    int id = 0;

    for(int idx = 1; idx < TStencil::Q; idx++) {

        leftorright = TStencil::Ci_xyz(0)[idx];

        if(leftorright == -1) {
        
            MPI_Isend(&obj.getDistribution()[mI0SendDistr[0]*TStencil::Q + idx],
                      1, distributionType, mNeighbors[0],
                      idx, MPI_COMM_WORLD, &commdist_request[id++]);
            MPI_Irecv(&obj.getDistribution()[mI0RecvDistr[1]*TStencil::Q + idx],
                      1, distributionType, mNeighbors[1],
                      idx, MPI_COMM_WORLD, &commdist_request[id++]);

        }
        else if(leftorright == 1) {
        
            MPI_Isend(&obj.getDistribution()[mI0SendDistr[1]*TStencil::Q + idx],
                      1, distributionType, mNeighbors[1],
                      idx, MPI_COMM_WORLD, &commdist_request[id++]);
            MPI_Irecv(&obj.getDistribution()[mI0RecvDistr[0]*TStencil::Q + idx],
                      1, distributionType, mNeighbors[0],
                      idx, MPI_COMM_WORLD, &commdist_request[id++]);

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

    using Stencil = typename TDistribution::Stencil;

    #pragma omp master
    {
    int nNeighbors = mNeighbors.size();
    MPI_Request commrequest[2*nNeighbors];

    for (int iNeighbor=0; iNeighbor<nNeighbors; iNeighbor++) {
        int tag = iNeighbor;
        MPI_Isend(&obj.mv_Distribution[mI0Send[iNeighbor]*Stencil::Q],
                  TNumNeighbors * TLattice::LY * TLattice::LZ * Stencil::Q,
                  mpi_get_type<double>(),
                  mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2*iNeighbor]);

        tag = (iNeighbor%2==0) ? iNeighbor+1 : iNeighbor-1;
        MPI_Irecv(&obj.mv_Distribution[mI0Recv[iNeighbor]*Stencil::Q],
                  TNumNeighbors * TLattice::LY * TLattice::LZ * Stencil::Q,
                  mpi_get_type<double>(),
                  mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2*iNeighbor+1]);
    }

    MPI_Waitall(2*nNeighbors, commrequest, MPI_STATUSES_IGNORE);
    }
    #endif

}

template<class TDerived, int TNumNeighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,TNumNeighbors>::communicateDistributionAllOld(TDistribution& obj) {

    #ifdef MPIPARALLEL
    if (mpi.size == 1) return;

    using Stencil = typename TDistribution::Stencil;

    #pragma omp master
    {
    int nNeighbors = mNeighbors.size();
    MPI_Request commrequest[2*nNeighbors];

    for (int iNeighbor=0; iNeighbor<nNeighbors; iNeighbor++) {
        int tag = iNeighbor;
        MPI_Isend(&obj.mv_OldDistribution[mI0Send[iNeighbor]*Stencil::Q],
                  TNumNeighbors * TLattice::LY * TLattice::LZ * Stencil::Q,
                  mpi_get_type<double>(),
                  mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2*iNeighbor]);

        tag = (iNeighbor%2==0) ? iNeighbor+1 : iNeighbor-1;
        MPI_Irecv(&obj.mv_OldDistribution[mI0Recv[iNeighbor]*Stencil::Q],
                  TNumNeighbors * TLattice::LY * TLattice::LZ * Stencil::Q,
                  mpi_get_type<double>(),
                  mNeighbors[iNeighbor], tag, MPI_COMM_WORLD, &commrequest[2*iNeighbor+1]);
    }

    MPI_Waitall(2*nNeighbors, commrequest, MPI_STATUSES_IGNORE);
    }
    #endif

}

/**
 * \brief NoParallel is a dummy implementation of Parallel that does not split the lattice or perform any communication.
 */
class NoParallel : public Parallel<NoParallel,0> {
    public:
        template<class TLattice>
        void init() {};

        #ifdef MPIPARALLEL
        template<class TLattice, class TStencil>
        static MPI_Datatype createDistributionType() { return MPI_DOUBLE; };
        #endif
};



/**
 * \brief X_Parallel contains functions and data for MPI parallelisation divided evenly in the X direction.
 * This class contains functions to communicate parameters for gradient calculations and to communicate
 * distributions for streaming.
 */
template<int TNumNeighbors=1>
class X_Parallel : public Parallel<X_Parallel<TNumNeighbors>,TNumNeighbors> {
    public:

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


#ifdef MPIPARALLEL
template<int TNumNeighbors>
template<class TLattice, class TStencil>
MPI_Datatype X_Parallel<TNumNeighbors>::createDistributionType() {
    MPI_Datatype distributionType;
    MPI_Type_vector(TLattice::LZ * TLattice::LY, 1, TStencil::Q, mpi_get_type<double>(), &distributionType);
    MPI_Type_commit(&distributionType);
    return distributionType;
}
#endif

template<int TNumNeighbors>
template<class TLattice>
void X_Parallel<TNumNeighbors>::init() {
    if (mpi.size <= 1) return;

    // Split lattice
    TLattice::HaloSize = TLattice::LY * TLattice::LZ * this->mMaxNeighbors;
    if (TLattice::LX % mpi.size == 0) {
        TLattice::LXdiv = (TLattice::LX / mpi.size + 2 * TNumNeighbors);
        TLattice::LXMPIOffset = (TLattice::LXdiv-2*TNumNeighbors)*mpi.rank;
    }
    else{
        throw std::runtime_error(std::string("Currently, the size of the domain in the x direction must be divisible by the number of mpi ranks."));
    }
    /*
    else if (mpi.rank<LX%mpi.size){
        LXdiv=((LX-LX%mpi.size)/mpi.size+1+2*TNumNeighbors);
    }
    else {
        LXdiv=((LX-LX%mpi.size)/mpi.size+2*TNumNeighbors);
    }
    */
    TLattice::N = TLattice::LXdiv * TLattice::LY * TLattice::LZ;

    // Define neighbor processors
    int leftNeighbor = mpi.rank - 1;
    if (leftNeighbor == -1) leftNeighbor = mpi.size - 1;
    int rightNeighbor = mpi.rank + 1;
    if (rightNeighbor == mpi.size) rightNeighbor = 0;
    this->mNeighbors = {leftNeighbor, rightNeighbor};

    // Create communication objects
    this->mI0Send = std::vector<int>(2);
    this->mI0Send[0] = TNumNeighbors * TLattice::LY * TLattice::LZ;
    this->mI0Send[1] = TLattice::N - (2*TNumNeighbors)*TLattice::LY*TLattice::LZ;
    this->mI0Recv = std::vector<int>(2);
    this->mI0Recv[0] = 0;
    this->mI0Recv[1] = TLattice::N - TNumNeighbors*TLattice::LY*TLattice::LZ;

    this->mI0SendDistr = std::vector<int>(2);
    this->mI0SendDistr[0] = (TNumNeighbors-1) * TLattice::LY * TLattice::LZ;
    this->mI0SendDistr[1] = TLattice::N - TNumNeighbors*TLattice::LY*TLattice::LZ;
    this->mI0RecvDistr = std::vector<int>(2);
    this->mI0RecvDistr[0] = (TNumNeighbors) * TLattice::LY * TLattice::LZ;
    this->mI0RecvDistr[1] = TLattice::N - (TNumNeighbors+1)*TLattice::LY*TLattice::LZ;

    // Create MPI buffer
    #ifdef MPIPARALLEL
    const int bufSize = (TLattice::LY * TLattice::LZ * TNumNeighbors * (5 + 2) * 2 * 2 + 1000) * sizeof(double);

    if(bufSize > MPIBUFFERSIZE) {

        if(MPIBUFFERSIZE != 0) MPI_Buffer_detach(MPIBUFFER, &MPIBUFFERSIZE);
        if(MPIBUFFERSIZE != 0) delete[] MPIBUFFER;
        MPIBUFFER = new char[bufSize];
        MPI_Buffer_attach(MPIBUFFER, bufSize);
        MPIBUFFERSIZE = bufSize;

    }
    #endif
}
