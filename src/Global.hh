#ifndef GLOBAL_HEADER
#define GLOBAL_HEADER

//Global.hh: Global parameters for the simulation

constexpr double DT=1.0; //Timestep
constexpr double CS2=1.0/3.0; //Speed of sound squared
constexpr int LX=10; //Size of lattice in x direction
constexpr int LY=10; //Size of lattice in y direction
constexpr int LZ=1; //Size of lattice in z direction
constexpr int TIMESTEPS=50; //Number of timesteps
constexpr int N=LX*LY*LZ; //Number of lattice points

#endif