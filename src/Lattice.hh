#pragma once
#include "Service.hh"
#include "Parallel.hh"
#include "Data.hh"

template<template<class, class, class> class data, template<class, class, int> class parallel, int lx, int ly, int lz = 1>
struct LatticeProperties{

    constexpr LatticeProperties<data, parallel, lx, ly, lz>& operator=(const LatticeProperties<data, parallel, lx, ly, lz>&) {}

    static constexpr int m_LX = lx;
    static constexpr int m_LY = ly;
    static constexpr int m_LZ = lz;

    #ifdef MPIPARALLEL

    constexpr LatticeProperties(double DT = 1.0) { m_DT=DT; }

    template<typename Stencil>
    using ParallelType = parallel<LatticeProperties<data, parallel, lx, ly, lz>,Stencil, 1>; //!<Chosen MPI parallelisation method when MPI enabled.
    static_assert(std::is_base_of<Parallel<LatticeProperties<data, parallel, lx, ly, lz>,1>,parallel<LatticeProperties<data, parallel, lx, ly, lz>, D2Q9,1>>::value,"ERROR: Chosen parallelisation method is not a parallelisation class.");
    
    static int m_LXdiv;
    static int m_LXMPIOffset;
    static int m_HaloSize;
    static int m_N;

    #else

    constexpr LatticeProperties(double DT = 1.0) { m_DT=DT; }

    template<typename Stencil>
    using ParallelType=No_Parallel<LatticeProperties<data, parallel, lx, ly, lz>, Stencil, 1>; //!<Default parallelisation when MPI disabled (Just serial).
    static_assert(std::is_base_of<Parallel<LatticeProperties<data, parallel, lx, ly, lz>,1>,parallel<LatticeProperties<data, parallel, lx, ly, lz>, D2Q9,1>>::value,"ERROR: Chosen parallelisation method is not a parallelisation class.");

    static constexpr int m_LXdiv = m_LX;
    static constexpr int m_HaloSize = 0;
    static constexpr int m_N = m_LX * m_LY * m_LZ;

    #endif

    template<typename Stencil>
    using DataType = data<LatticeProperties<data, parallel, lx, ly, lz>,Stencil, ParallelType<Stencil>>;//!<This will change the "DataType" implementation, which will govern the access of non-local data
    static_assert(std::is_base_of<Data_Base<LatticeProperties<data, parallel, lx, ly, lz>, D2Q9, ParallelType<D2Q9>>,data<LatticeProperties<data, parallel, lx, ly, lz>, D2Q9, ParallelType<D2Q9>>>::value,"ERROR: Chosen data method is not a data class.");

    static constexpr int m_NDIM = 3 - (lx <= 1 || ly <= 1 || lz <=1 );
    static double m_DT;
    
};

template<template<class, class, class> class data, template<class, class, int> class parallel, int NDIM>
struct LatticePropertiesRuntime {

    constexpr LatticePropertiesRuntime(int lx, int ly, double DT = 1.0) {
        m_DT = DT;
        m_LX = lx;
        m_LY = ly;
        m_LZ = 1;
        m_N = lx * ly;
    }
    constexpr LatticePropertiesRuntime(int lx, int ly, int lz, double DT = 1.0) { 
        m_DT = DT;
        m_LX = lx;
        m_LY = ly;
        m_LZ = lz;
        m_N = lx * ly * lz;   
    }

    constexpr LatticePropertiesRuntime<data, parallel, NDIM>& operator=(const LatticePropertiesRuntime<data, parallel, NDIM>&){

        return *this;

    }

    static const int m_LX;
    static const int m_LY;
    static const int m_LZ;
    static int m_LXdiv;
    static int m_LXMPIOffset;
    static int m_HaloSize;
    static int m_N;
    static constexpr int m_NDIM = NDIM;
    static double m_DT;

    #ifdef MPIPARALLEL

    template<typename Stencil>
    using ParallelType = parallel<LatticePropertiesRuntime<data, parallel, NDIM>,Stencil, 1>; //!<Chosen MPI parallelisation method when MPI enabled.
    static_assert(std::is_base_of<Parallel<LatticePropertiesRuntime<data, parallel, NDIM>,1>,parallel<LatticePropertiesRuntime<data, parallel, NDIM>,D2Q9,1>>::value,"ERROR: Chosen parallelisation method is not a parallelisation class.");

    #else

    template<typename Stencil>
    using ParallelType = No_Parallel<Stencil, 1>; //!<Default parallelisation when MPI disabled (Just serial).
    static_assert(std::is_base_of<Parallel<LatticePropertiesRuntime<data, parallel, NDIM>,1>,parallel<LatticePropertiesRuntime<data, parallel, NDIM>,D2Q9,1>>::value,"ERROR: Chosen parallelisation method is not a parallelisation class.");

    #endif

    template<typename Stencil>
    using DataType = data<LatticePropertiesRuntime<data, parallel, NDIM>,Stencil, ParallelType<Stencil>>;//!<This will change the "DataType" implementation, which will govern the access of non-local data
    static_assert(std::is_base_of<Data_Base<LatticePropertiesRuntime<data, parallel, NDIM>, D2Q9, ParallelType<D2Q9>>,data<LatticePropertiesRuntime<data, parallel, NDIM>, D2Q9, ParallelType<D2Q9>>>::value,"ERROR: Chosen data method is not a data class.");

};

template<template<class, class, class> class data, template<class, class, int> class parallel, int NDIM>
int LatticePropertiesRuntime<data, parallel, NDIM>::m_LXdiv;

template<template<class, class, class> class data, template<class, class, int> class parallel, int NDIM>
int LatticePropertiesRuntime<data, parallel, NDIM>::m_LXMPIOffset=0;

template<template<class, class, class> class data, template<class, class, int> class parallel, int NDIM>
int LatticePropertiesRuntime<data, parallel, NDIM>::m_HaloSize;

template<template<class, class, class> class data, template<class, class, int> class parallel, int NDIM>
int LatticePropertiesRuntime<data, parallel, NDIM>::m_N;

template<template<class, class, class> class data, template<class, class, int> class parallel, int NDIM>
double LatticePropertiesRuntime<data, parallel, NDIM>::m_DT=1.0;

template<template<class, class, class> class data, template<class, class, int> class parallel, int lx, int ly, int lz>
int LatticeProperties<data, parallel, lx, ly, lz>::m_LXdiv=lx;

template<template<class, class, class> class data, template<class, class, int> class parallel, int lx, int ly, int lz>
int LatticeProperties<data, parallel, lx, ly, lz>::m_LXMPIOffset=0;

template<template<class, class, class> class data, template<class, class, int> class parallel, int lx, int ly, int lz>
int LatticeProperties<data, parallel, lx, ly, lz>::m_HaloSize=0;

template<template<class, class, class> class data, template<class, class, int> class parallel, int lx, int ly, int lz>
int LatticeProperties<data, parallel, lx, ly, lz>::m_N=lx*ly*lz;

template<template<class, class, class> class data, template<class, class, int> class parallel, int lx, int ly, int lz>
double LatticeProperties<data, parallel, lx, ly, lz>::m_DT=1.0;

