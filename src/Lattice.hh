#pragma once
#include "Service.hh"
//template<int lx, int ly,int lz=1, template<class,template<class> class> class datatype=Data1, template<class,int> class paralleltype=X_Parallel>
template<template<class,class> class data,template<class,int> class parallel,int lx, int ly,int lz=1>
struct LatticeProperties{
    LatticeProperties(double DT=1.0):m_N(lx*ly*lz),m_DT(DT){
   
    }
    static constexpr int m_LX=lx;
    int m_LXdiv;
    int m_HaloSize;
    static constexpr int m_LY=ly;
    static constexpr int m_LZ=lz;
    int m_N;
    static constexpr int m_NDIM=(lx<=1||ly<=1||lz<=1)*-1+3;
    const double m_DT;

    #ifdef MPIPARALLEL
    template<typename Stencil>
    using ParallelType=parallel<Stencil,1>; //!< Chosen MPI parallelisation method when MPI enabled.
    #else
    template<typename... Stencil>
    using ParallelType=No_Parallel; //!< Default parallelisation when MPI disabled (Just serial).
    #endif
    template<typename Stencil>
    using DataType=data<Stencil,ParallelType<Stencil>>;//!< This will change the "DataType" implementation, which will govern the access of non-local data
};