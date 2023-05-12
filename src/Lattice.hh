#pragma once
#include "Service.hh"
#include "Parallel.hh"
#include "Data.hh"

template<template<class, class> class data, template<class, int> class parallel, int lx, int ly, int lz = 1>
struct LatticeProperties{

    #ifdef MPIPARALLEL

    constexpr LatticeProperties(double DT = 1.0) : m_N(lx * ly * lz), m_DT(DT) {}
    constexpr LatticeProperties(LatticeProperties<data, parallel, lx, ly, lz>& other) : m_N(other.m_N), m_DT(other.m_DT) {}

    #else

    constexpr LatticeProperties(double DT = 1.0):  m_DT(DT) {}
    constexpr LatticeProperties(LatticeProperties<data, parallel, lx, ly, lz>& other) : m_DT(other.m_DT) {}

    #endif

    constexpr LatticeProperties<data, parallel, lx, ly, lz>& operator=(const LatticeProperties<data, parallel, lx, ly, lz>&) {}

    static constexpr int m_LX = lx;
    static constexpr int m_LY = ly;
    static constexpr int m_LZ = lz;

    #ifdef MPIPARALLEL

    template<typename Stencil>
    using ParallelType = parallel<Stencil, 1>; //!<Chosen MPI parallelisation method when MPI enabled.
    //static_assert(std::is_base_of<Parallel<1>,parallel<D2Q9,1>>::value,"ERROR: Chosen parallelisation method is not a parallelisation class.");
    
    int m_LXdiv = m_LX;
    int m_HaloSize = 0;
    int m_N = lx * ly * lz;

    #else

    template<typename Stencil>
    using ParallelType=No_Parallel<Stencil, 1>; //!<Default parallelisation when MPI disabled (Just serial).

    static constexpr int m_LXdiv = m_LX;
    static constexpr int m_HaloSize = 0;
    static constexpr int m_N = m_LX * m_LY * m_LZ;

    #endif

    template<typename Stencil>
    using DataType = data<Stencil, ParallelType<Stencil>>;//!<This will change the "DataType" implementation, which will govern the access of non-local data
    //static_assert(std::is_base_of<Data_Base<D2Q9, ParallelType<D2Q9>>,data<D2Q9, ParallelType<D2Q9>>>::value,"ERROR: Chosen data method is not a data class.");

    static constexpr int m_NDIM = 3 - (lx <= 1 || ly <= 1 || lz <=1 );
    const double m_DT;
    
};

template<template<class, class> class data, template<class, int> class parallel, int NDIM>
struct LatticePropertiesRuntime {

    constexpr LatticePropertiesRuntime(int lx, int ly, double DT = 1.0) : m_LX(lx), m_LY(ly), m_LZ(1), m_N(lx*ly), m_DT(DT) {
   
    }
    constexpr LatticePropertiesRuntime(int lx, int ly, int lz, double DT = 1.0) : m_LX(lx), m_LY(ly), m_LZ(lz), m_N(lx * ly * lz), m_DT(DT) {
   
    }
    constexpr LatticePropertiesRuntime(LatticePropertiesRuntime<data, parallel, NDIM>& other) : m_LX(other.m_LX), m_LY(other.m_LY), m_LZ(other.m_LZ), m_N(other.m_N), m_DT(other.m_DT) {
   
    }
    constexpr LatticePropertiesRuntime<data, parallel, NDIM>& operator=(const LatticePropertiesRuntime<data, parallel, NDIM>&){

        return *this;

    }

    const int m_LX;
    const int m_LY;
    const int m_LZ;
    int m_LXdiv;
    int m_HaloSize;
    int m_N;
    static constexpr int m_NDIM = NDIM;
    const double m_DT;

    #ifdef MPIPARALLEL

    template<typename Stencil>
    using ParallelType = parallel<Stencil, 1>; //!<Chosen MPI parallelisation method when MPI enabled.
    //static_assert(std::is_base_of<Parallel<1>,parallel<D2Q9,1>>::value,"ERROR: Chosen parallelisation method is not a parallelisation class.");

    #else

    template<typename Stencil>
    using ParallelType = No_Parallel<Stencil, 1>; //!<Default parallelisation when MPI disabled (Just serial).

    #endif

    template<typename Stencil>
    using DataType = data<Stencil, ParallelType<Stencil>>;//!<This will change the "DataType" implementation, which will govern the access of non-local data
    //static_assert(std::is_base_of<Data_Base<D2Q9, ParallelType<D2Q9>>,data<D2Q9, ParallelType<D2Q9>>>::value,"ERROR: Chosen data method is not a data class.");

};
