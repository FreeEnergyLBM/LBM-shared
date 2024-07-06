#pragma once
#include <any>
#include <map>
#include <typeindex>
#include <typeinfo>
#include <utility>

#include "Data.hh"
#include "Parallel.hh"
#include "Service.hh"
#include "Stencil.hh"
#include "Template.hh"

/**
 * \file Lattice.hh
 * \brief This file contains classes which will hold information for the lattice.
 * Within these classes, we have fields for the lattice dimensions and information for the halo regions to be
 * communicated with MPI. These classes also contain functions to communicate parameters. Currently there are two
 * lattice classes, one in which the lattice dimensions are known at compile time and another where they will be passed
 * at runtime.
 */

/**
 * \brief This class will contain lattice information and handle parallelisation when the domain size is known at
 * compile time. \tparam TParallel Class to handle MPI parallelisation. \tparam lx Size of domain in x direction.
 * \tparam ly Size of domain in y direction. Defaults to 1.
 * \tparam lz Size of domain in z direction. Defaults to 1
 */
template <class TParallel, int lx, int ly = 1, int lz = 1>
struct LatticeProperties {
    using TLattice = LatticeProperties<TParallel, lx, ly, lz>;

    inline static void init() {}

    static constexpr int NDIM =
        3 - (lx <= 1 || ly <= 1 || lz <= 1) *
                (1 + ((lx <= 1 && ly <= 1) || (lx <= 1 && lz <= 1) ||
                      (ly <= 1 && lz <= 1)));              //!< Number of cartesian directions for the lattice.
    static constexpr int LX = lx;                          //!< Global length x.
    static constexpr int LY = ly;                          //!< Global length y.
    static constexpr int LZ = lz;                          //!< Global length z.
    inline static int N = lx * ly * lz;                    //!< Total size of local region (including halo).
    inline static int LXdiv = lx, LYdiv = ly, LZdiv = lz;  //!< Local lengths (including halo).
    inline static int LXMPIOffset = 0, LYMPIOffset = 0,
                      LZMPIOffset = 0;  //!< Global coordinates for the start of local regions (excluding halo).
    inline static int HaloXWidth = 0, HaloYWidth = 0, HaloZWidth = 0;  //!< Widths of the halo regions.
    inline static int HaloSize = 0;                //!< Size of each halo region (in each direction).
    inline static int subArray[3] = {lx, ly, lz};  //!< Local lengths (excluding halo).
    inline static double DT = 1;                   //!< Lattice imestep.
    inline static int Face[6] = {
        0, 0, 0, 0, 0, 0};        //!< Size of one slice parallel region. For 1D parallelisation only Face[0] is used.
    inline static int Width = 0;  //!< Chosen width of the parallel regions.
    static TParallel Parallel;    //!< Parallel class (contains communication functions).

    /**
     * \brief Initialises parallelisation and timestep variable.
     * \param DT Lattice timestep.
     */
    constexpr LatticeProperties(double DT = 1.0) {
        TLattice::DT = DT;
        Parallel.template init<TLattice>();  // Initialise the parallelisation
        initMPIBoundary<TLattice>();
    }

    /**
     * \brief Allows initialisation from another lattice.
     */
    constexpr TLattice& operator=(const TLattice&) { return *this; }

    static std::map<std::type_index, bool>
        alreadycommunicatedparameter;  //!< Map of parameter type to true if it has been communicated and false if not.

    enum { stream = 0, all = 1, allequilibrium = 2, allold = 3 };

    static std::unordered_map<std::pair<int, double*>, bool, pair_hash>
        alreadycommunicateddistribution;  //!< Map of pair containing distribution communication type and pointer to
                                          //!< distribution array to true if it has already been communicated and false
                                          //!< if not.

    /**
     * Resets all entries in alreadycommunicatedparameter and alreadycommunicateddistribution.
     */
    static void ResetParallelTracking() {
#pragma omp master
        {
            for (auto& [_, value] : alreadycommunicatedparameter) value = false;
            for (auto& [_, value] : alreadycommunicateddistribution) value = false;
        }
#pragma omp barrier
    }

    /**
     * \brief This function communicates a parameter adjacent to the parallel boundary into the halo regions of the
     * parameter in all neighboring lattices. \tparam TParameter type of the parameter. \param obj object of the
     * parameter.
     */
    template <class TParameter>
    static void communicate(TParameter& obj) {
#pragma omp master
        {
            // If parameter has not already been communicated or does not exist in the map.
            if (!alreadycommunicatedparameter.count(typeid(obj)) || !alreadycommunicatedparameter.at(typeid(obj))) {
                Parallel.template updateParameterBeforeCommunication<TLattice>(obj);
                Parallel.template communicateParameter<TLattice>(obj);
                Parallel.template updateParameterAfterCommunication<TLattice>(obj);
            }
            alreadycommunicatedparameter[typeid(obj)] = true;
        }
#pragma omp barrier
    }

    /**
     * \brief This function communicates multiple instances of a parameter
     * \tparam TParameter Type of object to be communicated.
     * \tparam NINS Number of instances to communicate
     * \tparam NDIR Number of directions in parameter
     */
    template <template <int> class TParameter, int NINS, int NDIR = 1>
    static void communicate();

    /**
     * \brief This function creates datatype for streaming distributions on its lattice.
     * \tparam TStencil The velocity stencil.
     */
    template <class TStencil>
    static void createDistributionType() {
        Parallel.template createDistributionType<TLattice, TStencil>();
    }

    /**
     * \brief This function destroys datatype for streaming distributions on its lattice.
     * \tparam TStencil The velocity stencil.
     */
    template <class TStencil>
    static void destroyDistributionType() {
        Parallel.template destroyDistributionType<TLattice>();
    }

    /**
     * \brief This function communicates the distributions that must be streamed across parallel boundaries.
     * \tparam TDistribution type of the distribution.
     * \param obj object of the distribution.
     */
    template <class TDistribution>
    static void communicateDistribution(TDistribution& obj) {
#pragma omp master
        {
            // If the distribution has not already been communicated with this communication type or does not exist in
            // the map.
            if (!alreadycommunicateddistribution.count(std::make_pair(0, obj.getDistributionPointer(0))) ||
                !alreadycommunicateddistribution.at(std::make_pair(0, obj.getDistributionPointer(0)))) {
                Parallel.template updateDistributionBeforeCommunication<TLattice>(obj);
                Parallel.template communicateDistribution<TLattice>(obj);
                Parallel.template updateDistributionAfterCommunication<TLattice>(obj);
            }
            alreadycommunicateddistribution[std::make_pair(0, obj.getDistributionPointer(0))] = true;
        }
#pragma omp barrier
    }

    /**
     * \brief This function communicates all distributions adjacent to the parallel boundary into the halo regions of
     * the distribution array in all neighboring lattices. \tparam TDistribution type of the distribution. \param obj
     * object of the distribution.
     */
    template <class TDistribution>
    static void communicateDistributionAll(TDistribution& obj) {
#pragma omp master
        {
            // If the distribution has not already been communicated with this communication type or does not exist in
            // the map.
            if (!alreadycommunicateddistribution.count(std::make_pair(1, obj.getDistributionPointer(0))) ||
                !alreadycommunicateddistribution.at(std::make_pair(1, obj.getDistributionPointer(0)))) {
                Parallel.template updateDistributionBeforeCommunicationAll<TLattice>(obj);
                Parallel.template communicateDistributionAll<TLattice>(obj);
                Parallel.template updateDistributionAfterCommunicationAll<TLattice>(obj);
            }
            alreadycommunicateddistribution[std::make_pair(1, obj.getDistributionPointer(0))] = true;
        }
#pragma omp barrier
    }

    /**
     * \brief This function communicates all equilibrium distributions adjacent to the parallel boundary into the halo
     * regions of the distribution array in all neighboring lattices. \tparam TDistribution type of the distribution.
     * \param obj object of the distribution.
     */
    template <class TDistribution>
    static void communicateDistributionAllEquilibrium(TDistribution& obj) {
#pragma omp master
        {
            // If the distribution has not already been communicated with this communication type or does not exist in
            // the map.
            if (!alreadycommunicateddistribution.count(std::make_pair(2, obj.getDistributionPointer(0))) ||
                !alreadycommunicateddistribution.at(std::make_pair(2, obj.getDistributionPointer(0)))) {
                Parallel.template updateDistributionBeforeCommunicationAllEquilibrium<TLattice>(obj);
                Parallel.template communicateDistributionAll<TLattice>(obj);
                Parallel.template updateDistributionAfterCommunicationAllEquilibrium<TLattice>(obj);
            }
            alreadycommunicateddistribution[std::make_pair(2, obj.getDistributionPointer(0))] = true;
        }
#pragma omp barrier
    }

    /**
     * \brief This function communicates all distributions from the previous timestep adjacent to the parallel boundary
     * into the halo regions of the distribution array in all neighboring lattices. \tparam TDistribution type of the
     * distribution. \param obj object of the distribution.
     */
    template <class TDistribution>
    static void communicateDistributionAllOld(TDistribution& obj) {
#pragma omp master
        {
            // If the distribution has not already been communicated with this communication type or does not exist in
            // the map.
            if (!alreadycommunicateddistribution.count(std::make_pair(3, obj.getDistributionPointer(0))) ||
                !alreadycommunicateddistribution.at(std::make_pair(3, obj.getDistributionPointer(0)))) {
                Parallel.template updateDistributionBeforeCommunicationAllOld<TLattice>(obj);
                Parallel.template communicateDistributionAll<TLattice>(obj);
                Parallel.template updateDistributionAfterCommunicationAllOld<TLattice>(obj);
            }
            alreadycommunicateddistribution[std::make_pair(3, obj.getDistributionPointer(0))] = true;
        }
#pragma omp barrier
    }
};

template <class TParallel, int lx, int ly, int lz>
TParallel LatticeProperties<TParallel, lx, ly, lz>::Parallel;

template <class TParallel, int lx, int ly, int lz>
std::map<std::type_index, bool> LatticeProperties<TParallel, lx, ly, lz>::alreadycommunicatedparameter;

template <class TParallel, int lx, int ly, int lz>
std::unordered_map<std::pair<int, double*>, bool, pair_hash>
    LatticeProperties<TParallel, lx, ly, lz>::alreadycommunicateddistribution;

template <class Lattice, class ParameterTuple, int NDIR>
struct communicateParameters;

template <class Lattice, class... Parameters, int NDIR>
struct communicateParameters<Lattice, std::tuple<Parameters...>, NDIR> {
    communicateParameters() { (Lattice::communicate(Parameters::template getInstance<Lattice, NDIR>()), ...); }
};

template <class TParallel, int lx, int ly, int lz>
template <template <int> class TParameter, int NINS, int NDIR>
void LatticeProperties<TParallel, lx, ly, lz>::communicate() {
    using Parameters = int_template<TParameter, int_sequence<NINS>>;
    communicateParameters<LatticeProperties<TParallel, lx, ly, lz>, Parameters, NDIR>();
}

/**
 * \brief This class will contain lattice information and handle parallelisation when the domain size is not known at
 * compile time. \tparam TParallel Class to handle MPI parallelisation. \tparam TNDIM Number of cartesian lattice
 * directions.
 */
template <class TParallel, int TNDIM>
struct LatticePropertiesRuntime {
    using TLattice = LatticePropertiesRuntime<TParallel, TNDIM>;

    //! Initialises parallel class.
    constexpr LatticePropertiesRuntime() { Parallel.template init<TLattice>(); }

    /**
     * \brief Initialises the lattice paramters with the given variables.
     * \param lx Size of domain in x direction.
     * \param ly Size of domain in y direction.
     * \param lz Size of domain in z direction.
     * \param DT Lattice timestep.
     */
    inline static void init(int lx, int ly, int lz, double DT = 1.0) {
        if (TNDIM <= 1 && ((lx > 1 && ly > 1) || (lx > 1 && lz > 1) || (ly > 1 && lz > 1)))
            throw std::runtime_error(
                "1D simulation but dimensions given to the LatticePropertiesRuntine class are not 1D");
        if (TNDIM <= 2 && (lx > 1 && ly > 1 && lz > 1))
            throw std::runtime_error(
                "2D simulation but dimensions given to the LatticePropertiesRuntine class are not 2D");
        LX = lx;
        LY = ly;
        LZ = lz;
        subArray[0] = lx;
        subArray[1] = ly;
        subArray[2] = lz;
        N = lx * ly * lz;
        LXdiv = lx;
        LYdiv = ly;
        LZdiv = lz;
        ;
        DT = DT;
    }

    //! Allows initialisation from another lattice.
    constexpr TLattice& operator=(const TLattice&) { return *this; }

    static constexpr int NDIM = TNDIM;                  //!< Number of cartesian lattice directions.
    inline static int LX = 0;                           //!< Global length x.
    inline static int LY = 0;                           //!< Global length y.
    inline static int LZ = 0;                           //!< Global length z.
    inline static int N = 0;                            //!< Total size of local region (including halo).
    inline static int LXdiv = 0, LYdiv = 0, LZdiv = 0;  //!< Local lengths (including halo).
    inline static int LXMPIOffset = 0, LYMPIOffset = 0,
                      LZMPIOffset = 0;  //!< Global coordinates for the start of local regions (excluding halo).
    inline static int HaloXWidth = 0, HaloYWidth = 0, HaloZWidth = 0;  //!< Widths of the halo regions.
    inline static int HaloSize = 0;      //!< Size of each halo region (in each direction).
    inline static int subArray[3] = {};  //!< Local lengths (excluding halo).
    inline static double DT = 1;         //!< Lattice imestep.
    inline static int Face[6] = {
        0, 0, 0, 0, 0, 0};        //!< Size of one slice parallel region. For 1D parallelisation only Face[0] is used.
    inline static int Width = 0;  //!< Chosen width of the parallel regions.
    static TParallel Parallel;    //!< Parallel class (contains communication functions).

    static std::map<std::type_index, bool>
        alreadycommunicatedparameter;  //!< Map of parameter type to true if it has been communicated and false if not.

    enum { stream = 0, all = 1, allequilibrium = 2, allold = 3 };

    static std::unordered_map<std::pair<int, double*>, bool, pair_hash>
        alreadycommunicateddistribution;  //!< Map of pair containing distribution communication type and pointer to
                                          //!< distribution array to true if it has already been communicated and false
                                          //!< if not.

    /**
     * Resets all entries in alreadycommunicatedparameter and alreadycommunicateddistribution.
     */
    static void ResetParallelTracking() {
#pragma omp master
        {
            for (auto& [_, value] : alreadycommunicatedparameter) value = false;
            for (auto& [_, value] : alreadycommunicateddistribution) value = false;
        }
#pragma omp barrier
    }

    /**
     * \brief This function communicates a parameter adjacent to the parallel boundary into the halo regions of the
     * parameter in all neighboring lattices. \tparam TParameter type of the parameter. \param obj object of the
     * parameter.
     */
    template <class TParameter>
    static void communicate(TParameter& obj) {
#pragma omp master
        {
            // If parameter has not already been communicated or does not exist in the map.
            if (!alreadycommunicatedparameter.count(typeid(obj)) || !alreadycommunicatedparameter.at(typeid(obj))) {
                Parallel.template updateParameterBeforeCommunication<TLattice>(obj);
                Parallel.template communicateParameter<TLattice>(obj);
                Parallel.template updateParameterAfterCommunication<TLattice>(obj);
            }
            alreadycommunicatedparameter[typeid(obj)] = true;
        }
#pragma omp barrier
    }

    //! This function creates datatype for streaming distributions on its lattice.
    template <class TStencil>
    static void createDistributionType() {
        Parallel.template createDistributionType<TLattice, TStencil>();
    }

    //! This function destroys datatype for streaming distributions on its lattice.
    template <class TStencil>
    static void destroyDistributionType() {
        Parallel.template destroyDistributionType<TLattice>();
    }

    /**
     * \brief This function communicates the distributions that must be streamed across parallel boundaries.
     * \tparam TDistribution type of the distribution.
     * \param obj object of the distribution.
     */
    template <class TDistribution>
    static void communicateDistribution(TDistribution& obj) {  // currently along X only

#pragma omp master
        {
            // If the distribution has not already been communicated with this communication type or does not exist in
            // the map.
            if (!alreadycommunicateddistribution.count(std::make_pair(0, obj.getDistributionPointer(0))) ||
                !alreadycommunicateddistribution.at(std::make_pair(0, obj.getDistributionPointer(0)))) {
                Parallel.template updateDistributionBeforeCommunication<TLattice>(obj);
                Parallel.template communicateDistribution<TLattice>(obj);
                Parallel.template updateDistributionAfterCommunication<TLattice>(obj);
            }
            alreadycommunicateddistribution[std::make_pair(0, obj.getDistributionPointer(0))] = true;
        }
#pragma omp barrier
    }

    /**
     * \brief This function communicates all distributions adjacent to the parallel boundary into the halo regions of
     * the distribution array in all neighboring lattices. \tparam TDistribution type of the distribution. \param obj
     * object of the distribution.
     */
    template <class TDistribution>
    static void communicateDistributionAll(TDistribution& obj) {
#pragma omp master
        {
            // If the distribution has not already been communicated with this communication type or does not exist in
            // the map.
            if (!alreadycommunicateddistribution.count(std::make_pair(1, obj.getDistributionPointer(0))) ||
                !alreadycommunicateddistribution.at(std::make_pair(1, obj.getDistributionPointer(0)))) {
                Parallel.template updateDistributionBeforeCommunicationAll<TLattice>(obj);
                Parallel.template communicateDistributionAll<TLattice>(obj);
                Parallel.template updateDistributionAfterCommunicationAll<TLattice>(obj);
            }
            alreadycommunicateddistribution[std::make_pair(1, obj.getDistributionPointer(0))] = true;
        }
#pragma omp barrier
    }

    /**
     * \brief This function communicates all eqiuilibrium distributions adjacent to the parallel boundary into the halo
     * regions of the distribution array in all neighboring lattices. \tparam TDistribution type of the distribution.
     * \param obj object of the distribution.
     */
    template <class TDistribution>
    static void communicateDistributionAllEquilibrium(TDistribution& obj) {
#pragma omp master
        {
            // If the distribution has not already been communicated with this communication type or does not exist in
            // the map.
            if (!alreadycommunicateddistribution.count(std::make_pair(2, obj.getDistributionPointer(0))) ||
                !alreadycommunicateddistribution.at(std::make_pair(2, obj.getDistributionPointer(0)))) {
                Parallel.template updateDistributionBeforeCommunicationAllEquilibrium<TLattice>(obj);
                Parallel.template communicateDistributionAll<TLattice>(obj);
                Parallel.template updateDistributionAfterCommunicationAllEquilibrium<TLattice>(obj);
            }
            alreadycommunicateddistribution[std::make_pair(2, obj.getDistributionPointer(0))] = true;
        }
#pragma omp barrier
    }

    /**
     * \brief This function communicates all distributions from the previous timestep adjacent to the parallel boundary
     * into the halo regions of the distribution array in all neighboring lattices. \tparam TDistribution type of the
     * distribution. \param obj object of the distribution.
     */
    template <class TDistribution>
    static void communicateDistributionAllOld(TDistribution& obj) {
#pragma omp master
        {
            // If the distribution has not already been communicated with this communication type or does not exist in
            // the map.
            if (!alreadycommunicateddistribution.count(std::make_pair(3, obj.getDistributionPointer(0))) ||
                !alreadycommunicateddistribution.at(std::make_pair(3, obj.getDistributionPointer(0)))) {
                Parallel.template updateDistributionBeforeCommunicationAllOld<TLattice>(obj);
                Parallel.template communicateDistributionAll<TLattice>(obj);
                Parallel.template updateDistributionAfterCommunicationAllOld<TLattice>(obj);
            }
            alreadycommunicateddistribution[std::make_pair(3, obj.getDistributionPointer(0))] = true;
        }
#pragma omp barrier
    }
};

template <class TParallel, int TNDIM>
TParallel LatticePropertiesRuntime<TParallel, TNDIM>::Parallel;

template <class TParallel, int TNDIM>
std::map<std::type_index, bool> LatticePropertiesRuntime<TParallel, TNDIM>::alreadycommunicatedparameter;

template <class TParallel, int TNDIM>
std::unordered_map<std::pair<int, double*>, bool, pair_hash>
    LatticePropertiesRuntime<TParallel, TNDIM>::alreadycommunicateddistribution;
