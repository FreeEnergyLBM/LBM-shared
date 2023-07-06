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
template<class TDerived, int num_neighbors=1>
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

    protected:
        int m_MaxNeighbors = 0;
        std::vector<int> m_Neighbors; //!<IDs of the neighboring processes.
        std::vector<int> m_I0Send; //!<Lattice index to begin sending to each processor.
        std::vector<int> m_I0Recv; //!<Lattice index to begin receiving from each processor.
        std::vector<int> m_I0SendDistr; //!<Lattice index to begin sending to each processor for the distributions.
        std::vector<int> m_I0RecvDistr; //!<Lattice index to begin receiving from each processor for the distributions.
};


template<class TDerived, int num_neighbors>
Parallel<TDerived,num_neighbors>::Parallel() {
    if(m_MaxNeighbors < num_neighbors) m_MaxNeighbors = num_neighbors;
}


/**
 * \details This will communicate the chosen parameter using MPI_Isend and MPI_Irecv, which are non-blocking methods of
 *          communication. This means that each process does not need to wait for the other processes to communicate. At
 *          the end of this function, we have a MPI_Waitall call, to ensure all processes are synced.
 */
template<class TDerived, int num_neighbors>
template<class TLattice, class TParameter>
inline void Parallel<TDerived,num_neighbors>::communicate(TParameter& obj) {
    #ifdef MPIPARALLEL
    if (mpi.size == 1) return;

    #pragma omp master
    {
    int nNeighbors = m_Neighbors.size();
    MPI_Request comm_request[2*nNeighbors];

    for (int iNeighbor=0; iNeighbor<nNeighbors; iNeighbor++) {
        int tag = iNeighbor;
        MPI_Isend(&obj.mv_Parameter[m_I0Send[iNeighbor]*obj.m_Num],
                  num_neighbors * TLattice::LY * TLattice::LZ * obj.m_Num,
                  mpi_get_type<typename TParameter::ParamType>(),
                  m_Neighbors[iNeighbor], tag, MPI_COMM_WORLD, &comm_request[2*iNeighbor]);

        tag = (iNeighbor%2==0) ? iNeighbor+1 : iNeighbor-1;
        MPI_Irecv(&obj.mv_Parameter[m_I0Recv[iNeighbor]*obj.m_Num],
                  num_neighbors * TLattice::LY * TLattice::LZ * obj.m_Num,
                  mpi_get_type<typename TParameter::ParamType>(),
                  m_Neighbors[iNeighbor], tag, MPI_COMM_WORLD, &comm_request[2*iNeighbor+1]);
    }

    MPI_Waitall(2*nNeighbors, comm_request, MPI_STATUSES_IGNORE);
    }
    #endif
}


template<class TDerived, int num_neighbors>
template<class TLattice, class TDistribution>
inline void Parallel<TDerived,num_neighbors>::communicateDistribution(TDistribution& obj) {

    #ifdef MPIPARALLEL
    if (mpi.size == 1) return;

    using Stencil = typename TDistribution::Stencil;
    static MPI_Datatype distributionType = TDerived::template createDistributionType<TLattice,Stencil>();

    #pragma omp master
    {
    MPI_Request comm_dist_request[20];

    int leftorright;
    int id = 0;

    for(int idx = 1; idx < Stencil::Q; idx++) {

        leftorright = Stencil::Ci_xyz(0)[idx];

        if(leftorright == -1) {
        
            MPI_Isend(&obj.getDistribution()[m_I0SendDistr[0]*Stencil::Q + idx],
                      1, distributionType, m_Neighbors[0],
                      idx, MPI_COMM_WORLD, &comm_dist_request[id++]);
            MPI_Irecv(&obj.getDistribution()[m_I0RecvDistr[1]*Stencil::Q + idx],
                      1, distributionType, m_Neighbors[1],
                      idx, MPI_COMM_WORLD, &comm_dist_request[id++]);

        }
        else if(leftorright == 1) {
        
            MPI_Isend(&obj.getDistribution()[m_I0SendDistr[1]*Stencil::Q + idx],
                      1, distributionType, m_Neighbors[1],
                      idx, MPI_COMM_WORLD, &comm_dist_request[id++]);
            MPI_Irecv(&obj.getDistribution()[m_I0RecvDistr[0]*Stencil::Q + idx],
                      1, distributionType, m_Neighbors[0],
                      idx, MPI_COMM_WORLD, &comm_dist_request[id++]);

        }

    }

    MPI_Waitall(id, comm_dist_request, MPI_STATUSES_IGNORE);
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
        template<class TLattice, class Stencil>
        static MPI_Datatype createDistributionType() { return MPI_DOUBLE; };
        #endif
};



/**
 * \brief X_Parallel contains functions and data for MPI parallelisation divided evenly in the X direction.
 * This class contains functions to communicate parameters for gradient calculations and to communicate
 * distributions for streaming.
 */
template<int num_neighbors=1>
class X_Parallel : public Parallel<X_Parallel<num_neighbors>,num_neighbors> {
    public:

        /**
         * \brief Initialise MPI variables for the parallelisation method.
         */
        template<class TLattice>
        void init();

        #ifdef MPIPARALLEL
        template<class TLattice, class Stencil>
        static MPI_Datatype createDistributionType(); //!<Create datatype for streaming distributions
        #endif
};


#ifdef MPIPARALLEL
template<int num_neighbors>
template<class TLattice, class Stencil>
MPI_Datatype X_Parallel<num_neighbors>::createDistributionType() {
    MPI_Datatype distributionType;
    MPI_Type_vector(TLattice::LZ * TLattice::LY, 1, Stencil::Q, mpi_get_type<double>(), &distributionType);
    MPI_Type_commit(&distributionType);
    return distributionType;
}
#endif

template<int num_neighbors>
template<class TLattice>
void X_Parallel<num_neighbors>::init() {
    if (mpi.size <= 1) return;

    // Split lattice
    TLattice::HaloSize = TLattice::LY * TLattice::LZ * this->m_MaxNeighbors;
    if (TLattice::LX % mpi.size == 0) {
        TLattice::LXdiv = (TLattice::LX / mpi.size + 2 * num_neighbors);
        TLattice::LXMPIOffset = (TLattice::LXdiv-2*num_neighbors)*mpi.rank;
    }
    else{
        throw std::runtime_error(std::string("Currently, the number of cores must be divisible by the size of the domain in the x direction."));
    }
    /*
    else if (mpi.rank<LX%mpi.size){
        LXdiv=((LX-LX%mpi.size)/mpi.size+1+2*num_neighbors);
    }
    else {
        LXdiv=((LX-LX%mpi.size)/mpi.size+2*num_neighbors);
    }
    */
    TLattice::N = TLattice::LXdiv * TLattice::LY * TLattice::LZ;

    // Define neighbor processors
    int leftNeighbor = mpi.rank - 1;
    if (leftNeighbor == -1) leftNeighbor = mpi.size - 1;
    int rightNeighbor = mpi.rank + 1;
    if (rightNeighbor == mpi.size) rightNeighbor = 0;
    this->m_Neighbors = {leftNeighbor, rightNeighbor};

    // Create communication objects
    this->m_I0Send = std::vector<int>(2);
    this->m_I0Send[0] = num_neighbors * TLattice::LY * TLattice::LZ;
    this->m_I0Send[1] = TLattice::N - (2*num_neighbors)*TLattice::LY*TLattice::LZ;
    this->m_I0Recv = std::vector<int>(2);
    this->m_I0Recv[0] = 0;
    this->m_I0Recv[1] = TLattice::N - num_neighbors*TLattice::LY*TLattice::LZ;

    this->m_I0SendDistr = std::vector<int>(2);
    this->m_I0SendDistr[0] = (num_neighbors-1) * TLattice::LY * TLattice::LZ;
    this->m_I0SendDistr[1] = TLattice::N - num_neighbors*TLattice::LY*TLattice::LZ;
    this->m_I0RecvDistr = std::vector<int>(2);
    this->m_I0RecvDistr[0] = (num_neighbors) * TLattice::LY * TLattice::LZ;
    this->m_I0RecvDistr[1] = TLattice::N - (num_neighbors+1)*TLattice::LY*TLattice::LZ;

    // Create MPI buffer
    #ifdef MPIPARALLEL
    const int bufSize = (TLattice::LY * TLattice::LZ * num_neighbors * (5 + 2) * 2 * 2 + 1000) * sizeof(double);

    if(bufSize > MPIBUFFERSIZE) {

        if(MPIBUFFERSIZE != 0) MPI_Buffer_detach(MPIBUFFER, &MPIBUFFERSIZE);
        if(MPIBUFFERSIZE != 0) delete[] MPIBUFFER;
        MPIBUFFER = new char[bufSize];
        MPI_Buffer_attach(MPIBUFFER, bufSize);
        MPIBUFFERSIZE = bufSize;

    }
    #endif
}
