#ifndef GLOBAL_HEADER
#define GLOBAL_HEADER
#ifdef MPIPARALLEL
#include <mpi.h>
#endif
//Global.hh: Global parameters for the simulation
#pragma omp begin declare target
constexpr int NO_NEIGHBOR=1;
constexpr double DT=1.0; //Timestep
constexpr double CS2=1.0/3.0; //Speed of sound squared
constexpr int LX=6400; //Size of lattice in x direction
constexpr int LY=3000; //Size of lattice in y direction
constexpr int LZ=1; //Size of lattice in z direction
constexpr int TIMESTEPS=10; //Number of timesteps
constexpr int SAVEINTERVAL=10;
#pragma omp end declare target
int N=LX*LY*LZ; //Number of lattice points
int LXdiv=LX;
int MAXNEIGHBORS=0;
int NUMPROCESSORS=1;
int CURPROCESSOR=0;
constexpr int NDIM=2;
#ifdef MPIPARALLEL
char* MPIBUFFER;
int MPIBUFFERSIZE;
MPI_Status status;
#endif
std::string DATA_DIR;
#ifdef OMPPARALLEL
double TOTALTIME=0;
#endif

#endif