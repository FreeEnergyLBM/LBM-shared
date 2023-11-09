#pragma once
#include "Service.hh"
#include "Parallel.hh"
#include "Data.hh"
#include <typeinfo>
#include <typeindex>
#include<map>

//====== LatticeProperties ======//

template<class TParallel, int lx, int ly, int lz = 1>
struct LatticeProperties{
    using TLattice = LatticeProperties<TParallel,lx,ly,lz>;

    inline static void init(){}

    //template<class TStencil>
    //using DataType = TData<TLattice, TStencil>;//!<This will change the "DataType" implementation, which will govern the access of non-local data
    //static_assert(std::is_base_of< Data_Base<TLattice,D1Q3>, DataType<D1Q3> >::value, "ERROR: Chosen data method is not a data class.");

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

    static std::map<std::type_index,bool> alreadycommunicatedparameter;

    enum {stream = 0, all = 1, allequilibrium = 2, allold = 3};

    static std::map<int,bool> alreadycommunicateddistribution;

    static void ResetParallelTracking() {
        #pragma omp master
        {
        for (auto& [_, value] : alreadycommunicatedparameter) value = false;
        for (auto& [_, value] : alreadycommunicateddistribution) value = false;
        }
    }

    //! This function communicates the halo regions of a parameter.
    template<class TParameter>
    static void communicate(TParameter& obj) {
        #pragma omp master
        {
        if (!alreadycommunicatedparameter[typeid(obj)]) {
            Parallel.template updateParameterBeforeCommunication<TLattice>(obj);
            Parallel.template communicateParameter<TLattice>(obj);
            Parallel.template updateParameterAfterCommunication<TLattice>(obj);
        }
        alreadycommunicatedparameter[typeid(obj)] = true;
        }
    }

    //! This function streams the distributions to the neighboring processor.
    template<class TDistribution>
    static void communicateDistribution(TDistribution& obj) { // currently along X only
        #pragma omp master
        {
        if (!alreadycommunicateddistribution[stream]) {
            Parallel.template updateDistributionBeforeCommunication<TLattice>(obj);
            Parallel.template communicateDistribution<TLattice>(obj);
            Parallel.template updateDistributionAfterCommunication<TLattice>(obj);
        }
        alreadycommunicateddistribution[stream] = true;
        }
    }

    template<class TDistribution>
    static void communicateDistributionAll(TDistribution& obj) {
        #pragma omp master
        {
        if (!alreadycommunicateddistribution[all]) {
            Parallel.template updateDistributionBeforeCommunicationAll<TLattice>(obj);
            Parallel.template communicateDistributionAll<TLattice>(obj);
            Parallel.template updateDistributionAfterCommunicationAll<TLattice>(obj);
        }
        alreadycommunicateddistribution[all] = true;
        }
    }

    template<class TDistribution>
    static void communicateDistributionAllEquilibrium(TDistribution& obj) {
        #pragma omp master
        {
        if (!alreadycommunicateddistribution[allequilibrium]) {
            Parallel.template updateDistributionBeforeCommunicationAllEquilibrium<TLattice>(obj);
            Parallel.template communicateDistributionAll<TLattice>(obj);
            Parallel.template updateDistributionAfterCommunicationAllEquilibrium<TLattice>(obj);
        }
        alreadycommunicateddistribution[allequilibrium] = true;
        }
    }

    template<class TDistribution>
    static void communicateDistributionAllOld(TDistribution& obj) {
        #pragma omp master
        {
        if (!alreadycommunicateddistribution[allold]) {
            Parallel.template updateDistributionBeforeCommunicationAllOld<TLattice>(obj);
            Parallel.template communicateDistributionAll<TLattice>(obj);
            Parallel.template updateDistributionAfterCommunicationAllOld<TLattice>(obj);
        }
        alreadycommunicateddistribution[allold] = true;
        }
    }
};

template<class TParallel, int lx, int ly, int lz>
TParallel LatticeProperties<TParallel,lx,ly,lz>::Parallel;

template<class TParallel, int lx, int ly, int lz>
std::map<std::type_index,bool> LatticeProperties<TParallel,lx,ly,lz>::alreadycommunicatedparameter;

template<class TParallel, int lx, int ly, int lz>
std::map<int,bool> LatticeProperties<TParallel,lx,ly,lz>::alreadycommunicateddistribution;

//====== LatticePropertiesRuntime ======//

template<class TParallel, int TNDIM>
struct LatticePropertiesRuntime {
    using TLattice = LatticePropertiesRuntime<TParallel,TNDIM>;

    constexpr LatticePropertiesRuntime() {
        Parallel.template init<TLattice>();
    }

    inline static void init(int lx, int ly, int lz, double DT=1.0){
        if(TNDIM<=1&&((lx>1&&ly>1)||(lx>1&&lz>1)||(ly>1&&lz>1))) throw std::runtime_error("1D simulation but dimensions given to the LatticePropertiesRuntine class are not 1D");
        if(TNDIM<=2&&(lx>1&&ly>1&&lz>1)) throw std::runtime_error("2D simulation but dimensions given to the LatticePropertiesRuntine class are not 2D");
        LX = lx;
        LY = ly;
        LZ = lz;
        subArray[0] = lx;
        subArray[1] = ly;
        subArray[2] = lz;
        N = lx * ly * lz;
        LXdiv = lx; LYdiv = ly; LZdiv = lz;;
        DT = DT;
    }

    constexpr TLattice& operator=(const TLattice&) {
        return *this;
    }

    static constexpr int NDIM = TNDIM;
    inline static int LX = 0;
    inline static int LY = 0;
    inline static int LZ = 0;
    inline static int N = 0;
    inline static int LXdiv = 0, LYdiv = 0, LZdiv = 0;
    inline static int LXMPIOffset = 0, LYMPIOffset = 0, LZMPIOffset = 0;
    inline static int HaloSize = 0, HaloXWidth = 0, HaloYWidth = 0, HaloZWidth = 0;
    inline static int subArray[3] = {};
    inline static double DT = 1.0;
    inline static int Face[6] = {}; // FaceX=0, FaceY=0, FaceZ=0, EdgeX=0, EdgeY=0, EdgeZ=0;
    inline static int Neighbors = 0; // FaceX=0, FaceY=0, FaceZ=0, EdgeX=0, EdgeY=0, EdgeZ=0;
    static TParallel Parallel;

    static std::map<std::type_index,bool> alreadycommunicatedparameter;

    enum {stream = 0, all = 1, allequilibrium = 2, allold = 3};

    static std::map<int,bool> alreadycommunicateddistribution;

    static void ResetParallelTracking() {
        #pragma omp master
        {
        for (auto& [_, value] : alreadycommunicatedparameter) value = false;
        for (auto& [_, value] : alreadycommunicateddistribution) value = false;
        }
    }

    //! This function communicates the halo regions of a parameter.
    template<class TParameter>
    static void communicate(TParameter& obj) {
        #pragma omp master
        {
        if (!alreadycommunicatedparameter[typeid(obj)]) {
            Parallel.template updateParameterBeforeCommunication<TLattice>(obj);
            Parallel.template communicateParameter<TLattice>(obj);
            Parallel.template updateParameterAfterCommunication<TLattice>(obj);
        }
        alreadycommunicatedparameter[typeid(obj)] = true;
        }
    }

    //! This function streams the distributions to the neighboring processor.
    template<class TDistribution>
    static void communicateDistribution(TDistribution& obj) { // currently along X only
        #pragma omp master
        {
        if (!alreadycommunicateddistribution[stream]) {
            Parallel.template updateDistributionBeforeCommunication<TLattice>(obj);
            Parallel.template communicateDistribution<TLattice>(obj);
            Parallel.template updateDistributionAfterCommunication<TLattice>(obj);
        }
        alreadycommunicateddistribution[stream] = true;
        }
    }

    template<class TDistribution>
    static void communicateDistributionAll(TDistribution& obj) {
        #pragma omp master
        {
        if (!alreadycommunicateddistribution[all]) {
            Parallel.template updateDistributionBeforeCommunicationAll<TLattice>(obj);
            Parallel.template communicateDistributionAll<TLattice>(obj);
            Parallel.template updateDistributionAfterCommunicationAll<TLattice>(obj);
        }
        alreadycommunicateddistribution[all] = true;
        }
    }

    template<class TDistribution>
    static void communicateDistributionAllEquilibrium(TDistribution& obj) {
        #pragma omp master
        {
        if (!alreadycommunicateddistribution[allequilibrium]) {
            Parallel.template updateDistributionBeforeCommunicationAllEquilibrium<TLattice>(obj);
            Parallel.template communicateDistributionAll<TLattice>(obj);
            Parallel.template updateDistributionAfterCommunicationAllEquilibrium<TLattice>(obj);
        }
        alreadycommunicateddistribution[allequilibrium] = true;
        }
    }

    template<class TDistribution>
    static void communicateDistributionAllOld(TDistribution& obj) {
        #pragma omp master
        {
        if (!alreadycommunicateddistribution[allold]) {
            Parallel.template updateDistributionBeforeCommunicationAllOld<TLattice>(obj);
            Parallel.template communicateDistributionAll<TLattice>(obj);
            Parallel.template updateDistributionAfterCommunicationAllOld<TLattice>(obj);
        }
        alreadycommunicateddistribution[allold] = true;
        }
    }
};

template<class TParallel, int TNDIM>
TParallel LatticePropertiesRuntime<TParallel,TNDIM>::Parallel;

template<class TParallel, int TNDIM>
std::map<std::type_index,bool> LatticePropertiesRuntime<TParallel,TNDIM>::alreadycommunicatedparameter;

template<class TParallel, int TNDIM>
std::map<int,bool> LatticePropertiesRuntime<TParallel,TNDIM>::alreadycommunicateddistribution;