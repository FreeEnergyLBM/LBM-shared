#pragma once
#ifdef MPIPARALLEL
#include <mpi.h>
#include "Global.hh"
#include "Service.hh"
#endif
/**
 * \file  Parallel.hh
 * \brief This contains classes to control the MPI parallelisation of the code.
 * The file contains a base class that will perform initialisation of the global MPI parameters. Other classes
 * inherit from this and provide funcions to communicate between processes.
 */

/**
 * \brief This class simply initialises global MPI parameters when constructed.
 * MAXNEIGHBORS is updated depending on the chosen number of neighbors. LXdiv (LX for each parallel block of
 * lattice points) is set based on the number of processors and number of neighbors chosen.
 */
template<int num_neighbors>
class Parallel {
    #ifdef MPIPARALLEL
    public:
        /**
         * \brief Constructor that updates global parameters. This is contained within a class as I may move stuff
         *        from X_Parallel to here.
         * MAXNEIGHBORS is updated depending on the chosen number of neighbors. LXdiv (LX for each parallel block of
         * lattice points) is set based on the number of processors and number of neighbors chosen. GETPROPERTIES().m_N is then
         * calculated as LXdiv*GETPROPERTIES().m_LY*GETPROPERTIES().m_LZ.
         */

        Parallel() {

            if(m_MaxNeighbors < num_neighbors) m_MaxNeighbors = num_neighbors;

            GETPROPERTIES().m_HaloSize = GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ * m_MaxNeighbors;

            if (GETPROPERTIES().m_LX % NUMPROCESSORS == 0) {
                GETPROPERTIES().m_LXdiv = (GETPROPERTIES().m_LX / NUMPROCESSORS + 2 * num_neighbors);
            }
            else{
                throw std::runtime_error(std::string("Currently, the number of cores must be divisible by the size of the domain in the x direction."));
            }
            /*
            else if (CURPROCESSOR<LX%NUMPROCESSORS){
                LXdiv=((LX-LX%NUMPROCESSORS)/NUMPROCESSORS+1+2*num_neighbors);
            }
            else {
                LXdiv=((LX-LX%NUMPROCESSORS)/NUMPROCESSORS+2*num_neighbors);
            }
            */
            GETPROPERTIES().m_N = GETPROPERTIES().m_LXdiv * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ;
            
        }

    private:

        int m_MaxNeighbors = 0;
    #endif
};

/**
 * \brief X_Parallel contains functions and data for MPI parallelisation divided evenly in the X direction.
 * This class contains functions to communicate parameters for gradient calculations and to communicate
 * distributions for streaming.
 */
template<class stencil, int num_neighbors>
class X_Parallel : public Parallel<num_neighbors> {
    #ifdef MPIPARALLEL
    public:

        /**
         * \brief Constructor that will initialise MPI variables for this parallelisation method.
         */

        X_Parallel();

        /**
         * \brief Function to fill halos of adjacent processors with the chosen parameter adjacent to the edge
         *        of the simulation in the X direction.
         * \param obj Object of chosen parameter.
         */
        template<class parameter>
        void communicate(parameter& obj);

        /**
         * \brief Function to update unknown distributions in the adjacent processors streamed from the edge of
         *        the current processor in the X direction.
         * \param obj Object of the distribution.
         */
        template<class distribution>
        void communicateDistribution(distribution& obj);

    private:

        int m_LeftNeighbor; //!<ID of the left neighbor of this process (in the X direction).
        int m_RightNeighbor; //!<ID of the right neighbor of this process (in the X direction).

        MPI_Datatype DistributionVector; //!<Datatype for streaming distributions (allows sending of one velocity index at a time) WILL NEED TO BE CHANGED BASED ON THE DATA TYPE.
    #endif
};
#ifdef MPIPARALLEL
/**
 * \details This will communicate the chosen parameter using MPI_Isend and MPI_Irecv, which are non-blocking methods of
 *          communication. This means that each process does not need to wait for the other processes to communicate. At
 *          the end of this function, we have a MPI_Waitall call, to ensure all processes are synced.
 */
template<class stencil, int num_neighbors>
template<class parameter>
void X_Parallel<stencil,num_neighbors>::communicate(parameter& obj) {

    MPI_Request comm_request[2];
    
    MPI_Isend(&obj.getParameter()[GETPROPERTIES().m_N * obj.m_Num - (num_neighbors + 1) * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ],
              num_neighbors * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ * obj.m_Num,
              mpi_get_type<typename parameter::ParamType>(),
              m_RightNeighbor, 0, MPI_COMM_WORLD, &comm_request[0]);

    MPI_Isend(&obj.getParameter()[num_neighbors * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ * obj.m_Num],
              num_neighbors * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ * obj.m_Num,
              mpi_get_type<typename parameter::ParamType>(),
              m_LeftNeighbor, 1, MPI_COMM_WORLD, &comm_request[1]);

    MPI_Irecv(&obj.getParameter()[GETPROPERTIES().m_N * obj.m_Num - num_neighbors * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ * obj.m_Num],
              num_neighbors * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ * obj.m_Num,
              mpi_get_type<typename parameter::ParamType>(),
              m_RightNeighbor, 1, MPI_COMM_WORLD, &comm_request[0]);

    MPI_Irecv(&obj.getParameter()[0],
              num_neighbors * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ * obj.m_Num,
              mpi_get_type<typename parameter::ParamType>(),
              m_LeftNeighbor, 0, MPI_COMM_WORLD, &comm_request[1]);
    
    MPI_Waitall(2, comm_request, MPI_STATUSES_IGNORE);

}

template<class stencil, int num_neighbors>
template<class distribution>
void X_Parallel<stencil,num_neighbors>::communicateDistribution(distribution& obj) {

    MPI_Request comm_dist_request[20];

    if (CURPROCESSOR < NUMPROCESSORS){

        int leftorright;
        int id = 0;

        for(int idx = 1; idx < stencil::Q; idx++) {

            leftorright = stencil::Ci_xyz(0)[idx];

            if(leftorright == -1) {
            
                MPI_Isend(&obj.getDistribution()[(num_neighbors-1) * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ * stencil::Q + idx],
                          1,
                          DistributionVector,
                          m_LeftNeighbor, id, MPI_COMM_WORLD, &comm_dist_request[id]);
                MPI_Irecv(&obj.getDistribution()[GETPROPERTIES().m_N * stencil::Q - (num_neighbors + 1) * GETPROPERTIES().m_LY*GETPROPERTIES().m_LZ * stencil::Q + idx],
                          1,
                          DistributionVector,
                          m_RightNeighbor, id, MPI_COMM_WORLD, &comm_dist_request[id]);
                id += 1;

            }
            else if(leftorright == 1) {
            
                MPI_Isend(&obj.getDistribution()[GETPROPERTIES().m_N * stencil::Q - (num_neighbors) * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ * stencil::Q + idx],
                          1,
                          DistributionVector,
                          m_LeftNeighbor, id, MPI_COMM_WORLD, &comm_dist_request[id]);
                MPI_Irecv(&obj.getDistribution()[(num_neighbors) * GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ * stencil::Q + idx],
                          1,
                          DistributionVector,
                          m_RightNeighbor, id, MPI_COMM_WORLD, &comm_dist_request[id]);
                id += 1;

            }
            
        }

        MPI_Waitall(id, comm_dist_request, MPI_STATUSES_IGNORE);

    }
    
}

template<class stencil, int num_neighbors>
X_Parallel<stencil, num_neighbors>::X_Parallel() {

    const int bufSize = (GETPROPERTIES().m_LY * GETPROPERTIES().m_LZ * num_neighbors * (5 + 2) * 2 * 2 + 1000) * sizeof(double);
    
    if(bufSize > MPIBUFFERSIZE) {

        if(MPIBUFFERSIZE != 0) MPI_Buffer_detach(MPIBUFFER, &MPIBUFFERSIZE);
        if(MPIBUFFERSIZE != 0) delete[] MPIBUFFER;
        MPIBUFFER = new char[bufSize];
        MPI_Buffer_attach(MPIBUFFER, bufSize);
        MPIBUFFERSIZE = bufSize;

    }
    
    m_LeftNeighbor = CURPROCESSOR-  1;
    if (m_LeftNeighbor == -1) m_LeftNeighbor = NUMPROCESSORS - 1;
    m_RightNeighbor = CURPROCESSOR + 1;
    if (m_RightNeighbor == NUMPROCESSORS) m_RightNeighbor = 0;

    MPI_Type_vector(GETPROPERTIES().m_LZ * GETPROPERTIES().m_LY, 1, stencil::Q, mpi_get_type<double>(), &DistributionVector);
    MPI_Type_commit(&DistributionVector);
    
}
#endif
