#pragma once
#include "Service.hh"
#include "Parallel.hh"
#include "Data.hh"
// TODO: Fix cyclic include of Data.hh



//====== LatticeProperties ======//

template<template<class, class> class TData, class TParallel, int lx, int ly, int lz = 1>
struct LatticeProperties{
    using TLattice = LatticeProperties<TData,TParallel,lx,ly,lz>;

    template<class TStencil>
    using DataType = TData<TLattice, TStencil>;//!<This will change the "DataType" implementation, which will govern the access of non-local data
    static_assert(std::is_base_of< Data_Base<TLattice,D1Q3>, DataType<D1Q3> >::value, "ERROR: Chosen data method is not a data class.");

    static constexpr int NDIM = 3 - (lx <= 1 || ly <= 1 || lz <=1 ) * (1 + ((lx <= 1 && ly <=1)||(lx <= 1 && lz <=1)||(ly <= 1 && lz <=1)));
    static constexpr int LX = lx;
    static constexpr int LY = ly;
    static constexpr int LZ = lz;
    inline static int N = lx * ly * lz;
    inline static int LXdiv = lx, LYdiv = ly, LZdiv = lz;
    inline static int LXMPIOffset = 0, LYMPIOffset = 0, LZMPIOffset = 0;
    inline static int HaloSize = 0, HaloXWidth = 0, HaloYWidth = 0, HaloZWidth = 0;
    inline static int subArray[3] = {lx, ly, lz};
    inline static double DT = 1;
    inline static int Face[6] = {0, 0, 0, 0, 0, 0}; // FaceX=0, FaceY=0, FaceZ=0, EdgeX=0, EdgeY=0, EdgeZ=0;
    inline static int Neighbors = 0; // FaceX=0, FaceY=0, FaceZ=0, EdgeX=0, EdgeY=0, EdgeZ=0;
    static TParallel Parallel;

    constexpr LatticeProperties(double DT=1.0) {
        DT = DT;
        Parallel.template init<TLattice>(); // Initialise the parallelisation
    }

    constexpr TLattice& operator=(const TLattice&) {
        return *this;
    }

    //! This function communicates the halo regions of a parameter.
    template<class TParameter>
    static void communicate(TParameter& obj) {
        Parallel.template updateParameterBeforeCommunication<TLattice>(obj);
        Parallel.template communicateParameter<TLattice>(obj);
        Parallel.template updateParameterAfterCommunication<TLattice>(obj);
    }

    //! This function streams the distributions to the neighboring processor.
    template<class TDistribution>
    static void communicateDistribution(TDistribution& obj) { // currently along X only
        Parallel.template updateDistributionBeforeCommunication<TLattice>(obj);
        Parallel.template communicateDistribution<TLattice>(obj);
        Parallel.template updateDistributionAfterCommunication<TLattice>(obj);
    }

    template<class TDistribution>
    static void communicateDistributionAll(TDistribution& obj) {
        Parallel.template updateDistributionBeforeCommunicationAll<TLattice>(obj);
        Parallel.template communicateDistributionAll<TLattice>(obj);
        Parallel.template updateDistributionAfterCommunicationAll<TLattice>(obj);
    }

    template<class TDistribution>
    static void communicateDistributionAllOld(TDistribution& obj) {
        Parallel.template updateDistributionBeforeCommunicationAllOld<TLattice>(obj);
        Parallel.template communicateDistributionAll<TLattice>(obj);
        Parallel.template updateDistributionAfterCommunicationAllOld<TLattice>(obj);
    }
};

template<template<class, class> class TData, class TParallel, int lx, int ly, int lz>
TParallel LatticeProperties<TData,TParallel,lx,ly,lz>::Parallel;

//====== LatticePropertiesRuntime ======//

template<template<class, class> class TData, class TParallel, int TNDIM>
struct LatticePropertiesRuntime {
    using TLattice = LatticePropertiesRuntime<TData,TParallel,TNDIM>;

    constexpr LatticePropertiesRuntime(int lx, int ly, double DT=1.0) : LatticePropertiesRuntime(lx, ly, 1, DT) {}

    constexpr LatticePropertiesRuntime(int lx, int ly, int lz, double DT=1.0) {
        if(TNDIM<=2&&lz>1) throw std::runtime_error("lz cannot be greater than 1 for a 2D simulation");
        LX = lx;
        LY = ly;
        LZ = lz;
        N = lx * ly * lz;
        LXdiv = lx; LYdiv = ly; LZdiv = lz;
        LXMPIOffset = 0; LYMPIOffset = 0; LZMPIOffset = 0;
        HaloSize = 0; HaloXWidth = 0; HaloYWidth = 0; HaloZWidth = 0;
        DT = DT;
    }

    constexpr TLattice& operator=(const TLattice&) {
        return *this;
    }

    static constexpr int NDIM = TNDIM;
    static int LX;
    static int LY;
    static int LZ;
    static int N;
    static int LXdiv, LYdiv, LZdiv;
    static int LXMPIOffset, LYMPIOffset, LZMPIOffset;
    static int HaloSize, HaloXWidth, HaloYWidth, HaloZWidth;
    static double DT;

    template<class TStencil>
    using DataType = TData<TLattice, TStencil>;//!<This will change the "DataType" implementation, which will govern the access of non-local data

    static_assert(std::is_base_of< Data_Base<TLattice,D1Q3>, DataType<D1Q3> >::value,
                  "ERROR: Chosen data method is not a data class.");
};
