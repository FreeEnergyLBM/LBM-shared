#pragma once
#include "Service.hh"
#include "Parallel.hh"
#include "Data.hh"



//====== LatticeProperties ======//

template<template<class, class> class TData, class TParallel, int lx, int ly, int lz = 1>
struct LatticeProperties{
    using TLattice = LatticeProperties<TData,TParallel,lx,ly,lz>;

    template<typename Stencil>
    using DataType = TData<TLattice, Stencil>;//!<This will change the "DataType" implementation, which will govern the access of non-local data
    static_assert(std::is_base_of< Data_Base<TLattice,D2Q9>, DataType<D2Q9> >::value, "ERROR: Chosen data method is not a data class.");

    static constexpr int m_NDIM = 3 - (lx <= 1 || ly <= 1 || lz <=1 );
    static constexpr int m_LX = lx;
    static constexpr int m_LY = ly;
    static constexpr int m_LZ = lz;
    inline static int m_N = lx * ly * lz;
    inline static int m_LXdiv = lx;
    inline static int m_LXMPIOffset = 0;
    inline static int m_HaloSize = 0;
    inline static double m_DT = 1;
    inline static TParallel m_Parallel;

    constexpr LatticeProperties(double DT=1.0) {
        m_DT = DT;
        m_Parallel.template init<TLattice>(); // Initialise the parallelisation
    }

    constexpr TLattice& operator=(const TLattice&) {
        return *this;
    }

    //! This function communicates the halo regions of a parameter.
    template<class TParameter>
    static void communicate(TParameter& obj) {
        m_Parallel.template communicate<TLattice>(obj);
    }

    //! This function streams the distributions to the neighboring processor.
    template<class TDistribution>
    static void communicateDistribution(TDistribution& obj) {
        m_Parallel.template communicateDistribution<TLattice>(obj);
    }
};


//====== LatticePropertiesRuntime ======//

template<template<class, class> class TData, class TParallel, int NDIM>
struct LatticePropertiesRuntime {
    using TLattice = LatticePropertiesRuntime<TData,TParallel,NDIM>;

    constexpr LatticePropertiesRuntime(int lx, int ly, double DT=1.0) : LatticePropertiesRuntime(lx, ly, 1, DT) {}

    constexpr LatticePropertiesRuntime(int lx, int ly, int lz, double DT=1.0) {
        m_LX = lx;
        m_LY = ly;
        m_LZ = lz;
        m_N = lx * ly * lz;
        m_LXdiv = lx;
        m_LXMPIOffset = 0;
        m_HaloSize = 0;
        m_DT = DT;
    }

    constexpr TLattice& operator=(const TLattice&) {
        return *this;
    }

    static constexpr int m_NDIM = NDIM;
    static int m_LX;
    static int m_LY;
    static int m_LZ;
    static int m_N;
    static int m_LXdiv;
    static int m_LXMPIOffset;
    static int m_HaloSize;
    static double m_DT;

    template<typename Stencil>
    using DataType = TData<TLattice, Stencil>;//!<This will change the "DataType" implementation, which will govern the access of non-local data

    static_assert(std::is_base_of< Data_Base<TLattice,D2Q9>, DataType<D2Q9> >::value,
                  "ERROR: Chosen data method is not a data class.");
};
