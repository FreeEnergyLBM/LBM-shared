#pragma once
#include "Service.hh"

template<template<class,class> class data,template<class,int> class parallel,int lx, int ly,int lz=1>
struct LatticeProperties{

    constexpr LatticeProperties(double DT=1.0):m_N(lx*ly*lz),m_DT(DT){
   
    }
    constexpr LatticeProperties(LatticeProperties<data,parallel,lx,ly,lz>& other):m_N(other.m_N),m_DT(other.m_DT){
   
    }
    constexpr LatticeProperties<data,parallel,lx,ly,lz>& operator=(const LatticeProperties<data,parallel,lx,ly,lz>&){
        return *this;
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

template<template<class,class> class data,template<class,int> class parallel,int NDIM>
struct LatticePropertiesRuntime{

    constexpr LatticePropertiesRuntime(int lx,int ly,double DT=1.0):m_LX(lx),m_LY(ly),m_LZ(1),m_N(lx*ly),m_DT(DT){
   
    }
    constexpr LatticePropertiesRuntime(int lx,int ly,int lz,double DT=1.0):m_LX(lx),m_LY(ly),m_LZ(lz),m_N(lx*ly*lz),m_DT(DT){
   
    }
    constexpr LatticePropertiesRuntime(LatticePropertiesRuntime<data,parallel,NDIM>& other):m_LX(other.m_LX),m_LY(other.m_LY),m_LZ(other.m_LZ),m_N(other.m_N),m_DT(other.m_DT){
   
    }
    constexpr LatticePropertiesRuntime<data,parallel,NDIM>& operator=(const LatticePropertiesRuntime<data,parallel,NDIM>&){
        return *this;
    }

    const int m_LX;
    const int m_LY;
    const int m_LZ;
    int m_LXdiv;
    int m_HaloSize;
    int m_N;
    static constexpr int m_NDIM=NDIM;
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
