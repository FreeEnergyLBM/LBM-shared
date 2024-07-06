#pragma once
#include <any>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Data.hh"
#include "Global.hh"
#include "Stencil.hh"

// Parameters.hh: This file details how macroscopic quantities are stored and interacted with.
// The Parameter class contains a static vector and functions to return the parameters stored.
// The reason for this static vector is to ensure that, if I use "Density" in multiple
// classes, the value will be consistent between classes. Note that, for every template configuration, a new
// static variable is created, so I pass each class to itself as a template to ensure the parameters are unique.

template <class TObj, class TLattice, typename T,
          int TNum = 1>  // obj template will guarantee a unique instance of the class with its own
                         // static vector. I pass the class to itself to guarantee this
                         //
                         // T determines the type stored in the vector, TNum
                         // determines the number of directions for each lattice point (eg you might
                         // want two directions corresponding to velocity in x and velocity in y at
                         // each point.
class Parameter {
   public:
    std::map<int, bool> mmInitialised;  // Change this to std::set

    static constexpr int mNum = TNum;
    static constexpr int mDirections = mNum;
    static constexpr const char *mName = TObj::mName;

    using ParamType = T;

    inline void save(int t);

    std::vector<T> mv_Parameter;      // Static vector (Does not change between objects of the class)
    std::vector<T> mv_CommParameter;  // Parameter vector reordered for convenience of communication

    template <class, class, int>
    friend class ParameterSingleton;

    inline std::vector<T> &getParameter() {  // Get a vector containing the parameters in the communication region
        return mv_Parameter;
    }

    inline std::vector<T> &getCommParameter() {  // Get a vector containing the parameters in the communication region
        return mv_CommParameter;
    }

   private:
    Parameter() {
        TLattice();                                    // Ensure the lattice variables have been initialised
        mv_Parameter.resize(TNum * TLattice::N, T());  // Resize to the desired size
        mv_CommParameter.resize(TNum * 4 * TLattice::HaloSize, T());  // Resize to the size of the communication region
    }

    Parameter(const Parameter &other) {}
};

template <class TObj, typename T = double, int TInstance = 0>
class ParameterSingleton {
   public:
    static constexpr int instance = TInstance;
    using ParamType = T;

    template <class TLattice, int TNum = 1>
    static inline Parameter<TObj, TLattice, T, TNum> &getInstance() {
        static Parameter<TObj, TLattice, T, TNum> instance;
        return instance;
    };

    template <class TLattice, int TNum = 1>
    static inline std::vector<T> &get() {  // Returns vector containing the parameter
        return getInstance<TLattice, TNum>().mv_Parameter;
    }

    // Returns pointer to parameter value at index idx
    template <class TLattice, int TNum = 1>
    static inline T *getAddress(const int idx) {
        return &getInstance<TLattice, TNum>().mv_Parameter[idx];
    }

    // Returns pointer to parameter value at position k and direction dir
    template <class TLattice, int TNum = 1>
    static inline T *getAddress(int k, int dir) {
        return &getInstance<TLattice, TNum>().mv_Parameter[k * TNum + dir];
    }

    template <class TLattice, int TNum = 1>
    static inline T &get(const int idx) {  // Returns const parameter at index idx

        // static Parameter<TObj,TLattice,T,TNum>& instance = getInstance<TLattice,TNum>();
        return getInstance<TLattice, TNum>().mv_Parameter[idx];
    }

    // Returns const parameter value at position k and direction dir
    template <class TLattice, int TNum = 1>
    static inline T &get(int k, int dir) {
        // MAKE MORE EFFICIENT!!!

        // static Parameter<TObj,TLattice,T,TNum>& instance = getInstance<TLattice,TNum>();
        return getInstance<TLattice, TNum>().mv_Parameter[k * TNum + dir];
    }

    template <class TLattice, int TNum = 1>
    static inline void initialise(const T val, const int k, const int dir = 0) {
#pragma omp critical
        {
            int i = k * TNum + dir;
            if (!getInstance<TLattice, TNum>().mmInitialised.count(i)) {
                getInstance<TLattice, TNum>().mv_Parameter[i] = val;
                getInstance<TLattice, TNum>().mmInitialised.insert({i, true});
            } else {
                getInstance<TLattice, TNum>().mmInitialised.erase(i);
            }
        }
    }

    template <class TLattice, int TNum = 1, int dir = 0>
    static void smooth(int iter = 1) {
#pragma omp master
        {
            using Stencil =
                std::conditional_t<TLattice::NDIM == 1, D1Q3, std::conditional_t<TLattice::NDIM == 2, D2Q9, D3Q27>>;

            using DataType = Data_Base<TLattice, Stencil>;
            std::vector<T> vals;
            vals.resize(TNum * TLattice::N);
            DataType &data = DataType::getInstance();

            for (int i = 0; i < iter; i++) {
                TLattice::ResetParallelTracking();
                TLattice::communicate(getInstance<TLattice, TNum>());
                for (int k : RangeK<TLattice>()) {
                    int iData = k * TNum + dir;
                    getInstance<TLattice, TNum>().mmInitialised.erase(iData);
                    T sum = 0;
                    for (int idx = 0; idx < Stencil::Q; idx++) {
                        sum += get<TLattice, TNum>(data.getNeighbor(k, idx), dir);
                    }
                    vals[iData] = sum / Stencil::Q;
                }
                for (int k : RangeK<TLattice>()) {
                    initialise<TLattice, TNum>(vals[k * TNum + dir], k, dir);
                }
            }
            TLattice::ResetParallelTracking();
            TLattice::communicate(getInstance<TLattice, TNum>());
            TLattice::ResetParallelTracking();
        }
    }

    template <class TLattice, int TNum = 1, int dir = 0>
    static void set(T val) {
        for (int k : RangeK<TLattice>()) {
            initialise<TLattice, TNum>(val, k, dir);
        }
    }

    template <class TLattice, int TNum = 1, int dir = 0>
    static void set(T (*condition)(const int)) {
        for (int k : RangeK<TLattice>()) {
            initialise<TLattice, TNum>(condition(k), k, dir);
        }
    }

    template <class TLattice, int TNum = 1, int dir = 0>
    static void set(bool (*condition)(const int), T val) {
        for (int k : RangeK<TLattice>()) {
            if (condition(k)) initialise<TLattice, TNum>(val, k, dir);
        }
    }

    template <class TLattice, int TNum = 1, int dir = 0>
    static void set(bool (*condition)(const int), T val, T false_val) {
        for (int k : RangeK<TLattice>()) {
            if (condition(k)) {
                initialise<TLattice, TNum>(val, k, dir);
            } else {
                initialise<TLattice, TNum>(false_val, k, dir);
            }
        }
    }

    ParameterSingleton(ParameterSingleton<TObj, T, TInstance> const &) = delete;

    void operator=(ParameterSingleton<TObj, T, TInstance> const &) = delete;

   private:
    ParameterSingleton(){};
};

// Utility function to get a parameter instance at runtime
// e.g. getInstance<OrderParameter, 3, TLattice>(1) will return OrderParameter<1>::get<TLattice>()
// Recursion is needed to go from a run-time integer to a compile-time integer
template <template <int> typename TParameter, typename TLattice, int NDIM, typename T, T I, T... Is>
std::vector<double> &getInstance(int instance, std::index_sequence<I, Is...>) {
    if (instance == I) {
        return TParameter<I>::template get<TLattice, NDIM>();
    } else if constexpr (sizeof...(Is) > 0) {
        return getInstance<TParameter, TLattice, NDIM>(instance, std::index_sequence<Is...>{});
    } else {
        throw std::invalid_argument("Requested instance not in index list.");
    }
}

template <template <int> typename TParameter, int N, typename TLattice, int NDIM = 1>
std::vector<double> &getInstance(int instance) {
    if (instance >= N) throw std::invalid_argument("Requested instance must be less than N.");
    return getInstance<TParameter, TLattice, NDIM>(instance, std::make_index_sequence<N>{});
}

// Utility function to get a parameter gradient instance at runtime
// e.g. getGradientInstance<Gradient, OrderParameter, 3, TLattice>(1)
// will return Gradient<OrderParameter<1>>::get<TLattice>()
// Recursion is needed to go from a run-time integer to a compile-time integer
template <template <typename> typename TGradient, template <int> typename TParameter, typename TLattice, int NDIM,
          typename T, T I, T... Is>
std::vector<double> &getGradientInstance(int instance, std::index_sequence<I, Is...>) {
    if (instance == I) {
        return TGradient<TParameter<I>>::template get<TLattice, NDIM>();
    } else if constexpr (sizeof...(Is) > 0) {
        return getGradientInstance<TGradient, TParameter, TLattice, NDIM>(instance, std::index_sequence<Is...>{});
    } else {
        throw std::invalid_argument("Requested instance not in index list.");
    }
}

template <template <typename> typename TGradient, template <int> typename TParameter, int N, typename TLattice,
          int NDIM = 1>
std::vector<double> &getGradientInstance(int instance) {
    if (instance >= N) throw std::invalid_argument("Requested instance must be less than N.");
    return getGradientInstance<TGradient, TParameter, TLattice, NDIM>(instance, std::make_index_sequence<N>{});
}

template <class TLattice>
class SaveHandler;

template <class TObj, class TLattice, typename T, int TNum>
inline void Parameter<TObj, TLattice, T, TNum>::save(int t) {  // Function to save parameter stored in this class
    SaveHandler<TLattice>::template saveParameter<TObj, TNum>(t);
}

template <int TNDIM>
struct Boundary {
    int Id;
    bool IsCorner;
    std::array<int8_t, TNDIM> NormalDirection;
};

template <int TNDIM, int TInstance = 0>
struct BoundaryLabels : public ParameterSingleton<BoundaryLabels<TNDIM, TInstance>, Boundary<TNDIM>, TInstance> {
    static constexpr const char *mName = "BoundaryLabels";
};  // Labelling of geometry

template <int TInstance = 0>
struct Velocity : public ParameterSingleton<Velocity<TInstance>, double, TInstance> {
    static constexpr const char *mName = "Velocity";

};  // Velocity, with directions D corresponding to the number of cartesian directions in the stencil

template <int TInstance = 0>
struct VelocityOld : public ParameterSingleton<VelocityOld<TInstance>, double, TInstance> {
    static constexpr const char *mName = "VelocityOld";
};

template <int TInstance = 0>
struct ViscousStress : public ParameterSingleton<ViscousStress<TInstance>, double, TInstance> {
    static constexpr const char *mName = "ViscousStress";
};

template <int TInstance = 0>
struct ViscousDissipation : public ParameterSingleton<ViscousDissipation<TInstance>, double, TInstance> {
    static constexpr const char *mName = "ViscousDissipation";
};

template <int TInstance = 0>
struct Density : public ParameterSingleton<Density<TInstance>, double, TInstance> {
    static constexpr const char *mName = "Density";

};  // Density

template <int TInstance = 0>
struct DensityOld : public ParameterSingleton<DensityOld<TInstance>, double, TInstance> {
    static constexpr const char *mName = "DensityOld";

};  // Density

template <int TInstance = 0>
struct Pressure : public ParameterSingleton<Pressure<TInstance>, double, TInstance> {
    static constexpr const char *mName = "Pressure";

};  // Presure

template <int TInstance = 0>
struct PressureOld : public ParameterSingleton<PressureOld<TInstance>, double, TInstance> {
    static constexpr const char *mName = "PressureOld";

};  // Presure

template <int TInstance = 0>
struct OrderParameter : public ParameterSingleton<OrderParameter<TInstance>, double, TInstance> {
    static constexpr const char *mName = "OrderParameter";

};  // Order parameter representing relative concentration of the phases

template <int TInstance = 0>
struct ChemicalPotential : public ParameterSingleton<ChemicalPotential<TInstance>, double, TInstance> {
    static constexpr const char *mName = "ChemicalPotential";

};  // Chemical potential for the multicomponent model

template <int TInstance = 0>
struct ChemicalPotentialOld : public ParameterSingleton<ChemicalPotentialOld<TInstance>, double, TInstance> {
    static constexpr const char *mName = "ChemicalPotentialOld";

};  // Chemical potential for the multicomponent model

template <int TInstance = 0>
struct OrderParameterOld : public ParameterSingleton<OrderParameterOld<TInstance>, double, TInstance> {
    static constexpr const char *mName = "OrderParameterOld";
};

template <int TInstance = 0>
struct Humidity : public ParameterSingleton<Humidity<TInstance>> {
    static constexpr const char *mName = "Humidity";
};

template <int TInstance = 0>
struct Temperature : public ParameterSingleton<Temperature<TInstance>> {
    static constexpr const char *mName = "Humidity";
};

template <int TInstance = 0>
struct HumidityOld : public ParameterSingleton<HumidityOld<TInstance>> {
    static constexpr const char *mName = "HumidityOld";
};

template <int TInstance = 0>
struct MassSink : public ParameterSingleton<MassSink<TInstance>> {
    static constexpr const char *mName = "MassSink";
};

template <int TInstance = 0>
struct Solute : public ParameterSingleton<Solute<TInstance>> {
    static constexpr const char *mName = "Solute";
};

template <int TInstance = 0>
struct SoluteOld : public ParameterSingleton<SoluteOld<TInstance>> {
    static constexpr const char *mName = "SoluteOld";
};

template <class TObj>
struct Laplacian : public ParameterSingleton<Laplacian<TObj>, double> {
    static constexpr char mName[9 + sizeof(TObj::mName)] = "Laplacian" + TObj::mName;
};  // Directional first order gradients of the order parameter

template <int TInstance = 0>
struct LaplacianChemicalPotential : public Laplacian<ChemicalPotential<TInstance>> {
    static constexpr const char *mName = "LaplacianChemicalPotential";

};  // Laplacian of the order parameter

template <int TInstance = 0>
struct LaplacianDensity : public Laplacian<Density<TInstance>> {
    static constexpr const char *mName = "LaplacianDensity";

};  // Laplacian of the order parameter

template <int TInstance = 0>
struct LaplacianOrderParameter : public Laplacian<OrderParameter<TInstance>> {
    static constexpr const char *mName = "LaplacianOrderParameter";

};  // Laplacian of the order parameter

template <class TObj>
struct Gradient : public ParameterSingleton<Gradient<TObj>, double> {
    static constexpr char mName[8 + sizeof(TObj::mName)] = "Gradient" + TObj::mName;
    using ParameterSingleton<Gradient<TObj>, double>::get;
    template <class TLattice, int TNumGrad, int TNum>
    static inline typename TObj::ParamType &get(int k, int dir1, int dir2) {
        return ParameterSingleton<Gradient<TObj>, double>::template getInstance<TLattice, TNum * TNumGrad>()
            .mv_Parameter[k * TNumGrad * TNum + TNumGrad * dir1 + dir2];
    }
};  // Directional first order gradients of the order parameter

template <class TObj>
struct GradientMixed : public ParameterSingleton<GradientMixed<TObj>, double> {
    static constexpr char mName[13 + sizeof(TObj::mName)] = "GradientMixed" + TObj::mName;
    using ParameterSingleton<GradientMixed<TObj>, double>::get;
    template <class TLattice, int TNumGrad, int TNum>
    static inline typename TObj::ParamType &get(int k, int dir1, int dir2) {
        return ParameterSingleton<GradientMixed<TObj>, double>::template getInstance<TLattice, TNum * TNumGrad>()
            .mv_Parameter[k * TNumGrad * TNum + TNumGrad * dir1 + dir2];
    }
};  // Directional first order gradients of the order parameter

template <int TInstance = 0>
struct GradientVelocity : public Gradient<Velocity<TInstance>> {
    static constexpr const char *mName = "GradientVelocity";
};

template <int TInstance = 0>
struct GradientOrderParameter : public Gradient<OrderParameter<TInstance>> {
    static constexpr const char *mName = "GradientOrderParameter";

};  // Directional first order gradients of the order parameter

template <int TInstance = 0>
struct GradientDensity : public Gradient<Density<TInstance>> {
    static constexpr const char *mName = "GradientDensity";

};  // Directional first order gradients of the order parameter

template <int TInstance = 0>
struct GradientChemicalPotential : public Gradient<ChemicalPotential<TInstance>> {
    static constexpr const char *mName = "GradientChemicalPotential";

};  // Directional first order gradients of the order parameter

template <int TInstance = 0>
struct GradientPressure : public Gradient<Pressure<TInstance>> {
    static constexpr const char *mName = "GradientPressure";

};  // Directional first order gradients of the order parameter

template <int TInstance = 0>
struct GradientHumidity : public Gradient<Humidity<TInstance>> {
    static constexpr const char *mName = "GradientHumidity";
};

template <int TInstance = 0>
struct MixedGradientOrderParameter : public GradientMixed<OrderParameter<TInstance>> {
    static constexpr const char *mName = "MixedGradientOrderParameter";

};  // Directional first order gradients of the order parameter

template <int TInstance = 0>
struct MixedGradientDensity : public GradientMixed<Density<TInstance>> {
    static constexpr const char *mName = "MixedGradientDensity";

};  // Directional first order gradients of the order parameter

template <int TInstance = 0>
struct MixedGradientPressure : public GradientMixed<Pressure<TInstance>> {
    static constexpr const char *mName = "MixedGradientPressure";

};  // Directional first order gradients of the order parameter

template <int TInstance = 0>
struct Tau : public ParameterSingleton<Tau<TInstance>, double, TInstance> {
    static constexpr const char *mName = "Tau";
};  // Labelling of geometry

template <int TInstance = 0>
struct InverseTau : public ParameterSingleton<InverseTau<TInstance>, double, TInstance> {
    static constexpr const char *mName = "InverseTau";
};  // Labelling of geometry

template <int TInstance = 0>
struct LogFugacity : public ParameterSingleton<LogFugacity<TInstance>, double, TInstance> {
    static constexpr const char *mName = "Fugacity";
};  // Labelling of geometry

template <int TInstance = 0>
struct GradientLogFugacity : public ParameterSingleton<GradientLogFugacity<TInstance>, double, TInstance> {
    static constexpr const char *mName = "Fugacity";
};  // Labelling of geometry

#ifdef MPIPARALLEL
MPI_Datatype mMPIBoundary;
bool mMPIBoundaryInitialised = false;
#endif

template <class TLattice>
void initMPIBoundary() {
#ifdef MPIPARALLEL
    if (mMPIBoundaryInitialised) return;

    MPI_Type_create_resized(MPI_INT, 0L, sizeof(Boundary<TLattice::NDIM>), &mMPIBoundary);
    MPI_Type_commit(&mMPIBoundary);

    int blocklengths[3] = {1, 1, TLattice::NDIM};
    MPI_Datatype types[3] = {MPI_INT, MPI_CXX_BOOL, MPI_INT8_T};
    std::array<int8_t, TLattice::NDIM> a;
    MPI_Aint offsets[3];
    offsets[0] = offsetof(Boundary<TLattice::NDIM>, Id);
    offsets[1] = offsetof(Boundary<TLattice::NDIM>, IsCorner);
    offsets[2] =
        offsetof(Boundary<TLattice::NDIM>, NormalDirection) +
        (size_t)((((char *)(&a) - (char *)(&a[0]))));  // + offsetof(std::array<int8_t,Lattice::NDIM>, NormalDirection);

    MPI_Type_create_struct(3, blocklengths, offsets, types, &mMPIBoundary);
    MPI_Type_commit(&mMPIBoundary);
    mMPIBoundaryInitialised = true;
#endif
}

#ifdef MPIPARALLEL
/**\fn      mpi_get_type
 * \brief   Small template function to return the correct MPI_DATATYPE
 *          data type need for an MPI message as a constexpr at compile time
 *          https://www.mpich.org/static/docs/latest/www3/Constants.html
 *          Call in a template function with mpi_get_type<T>()
 *
 * \tparam  T   The C++ data type used in the MPI function
 * \return  The MPI_Datatype belonging to the template C++ data type T
 */
template <typename T, class TLattice>
[[nodiscard]] constexpr MPI_Datatype mpi_get_type() noexcept {
    MPI_Datatype mpi_type = MPI_DATATYPE_NULL;

    if constexpr (std::is_same_v<T, char>) {
        mpi_type = MPI_CHAR;
    } else if constexpr (std::is_same_v<T, signed char>) {
        mpi_type = MPI_SIGNED_CHAR;
    } else if constexpr (std::is_same_v<T, unsigned char>) {
        mpi_type = MPI_UNSIGNED_CHAR;
    } else if constexpr (std::is_same_v<T, wchar_t>) {
        mpi_type = MPI_WCHAR;
    } else if constexpr (std::is_same_v<T, signed short>) {
        mpi_type = MPI_SHORT;
    } else if constexpr (std::is_same_v<T, unsigned short>) {
        mpi_type = MPI_UNSIGNED_SHORT;
    } else if constexpr (std::is_same_v<T, signed int>) {
        mpi_type = MPI_INT;
    } else if constexpr (std::is_same_v<T, unsigned int>) {
        mpi_type = MPI_UNSIGNED;
    } else if constexpr (std::is_same_v<T, signed long int>) {
        mpi_type = MPI_LONG;
    } else if constexpr (std::is_same_v<T, unsigned long int>) {
        mpi_type = MPI_UNSIGNED_LONG;
    } else if constexpr (std::is_same_v<T, signed long long int>) {
        mpi_type = MPI_LONG_LONG;
    } else if constexpr (std::is_same_v<T, unsigned long long int>) {
        mpi_type = MPI_UNSIGNED_LONG_LONG;
    } else if constexpr (std::is_same_v<T, float>) {
        mpi_type = MPI_FLOAT;
    } else if constexpr (std::is_same_v<T, double>) {
        mpi_type = MPI_DOUBLE;
    } else if constexpr (std::is_same_v<T, long double>) {
        mpi_type = MPI_LONG_DOUBLE;
    } else if constexpr (std::is_same_v<T, std::int8_t>) {
        mpi_type = MPI_INT8_T;
    } else if constexpr (std::is_same_v<T, std::int16_t>) {
        mpi_type = MPI_INT16_T;
    } else if constexpr (std::is_same_v<T, std::int32_t>) {
        mpi_type = MPI_INT32_T;
    } else if constexpr (std::is_same_v<T, std::int64_t>) {
        mpi_type = MPI_INT64_T;
    } else if constexpr (std::is_same_v<T, std::uint8_t>) {
        mpi_type = MPI_UINT8_T;
    } else if constexpr (std::is_same_v<T, std::uint16_t>) {
        mpi_type = MPI_UINT16_T;
    } else if constexpr (std::is_same_v<T, std::uint32_t>) {
        mpi_type = MPI_UINT32_T;
    } else if constexpr (std::is_same_v<T, std::uint64_t>) {
        mpi_type = MPI_UINT64_T;
    } else if constexpr (std::is_same_v<T, bool>) {
        mpi_type = MPI_C_BOOL;
    } else if constexpr (std::is_same_v<T, std::complex<float>>) {
        mpi_type = MPI_C_COMPLEX;
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        mpi_type = MPI_C_DOUBLE_COMPLEX;
    } else if constexpr (std::is_same_v<T, std::complex<long double>>) {
        mpi_type = MPI_C_LONG_DOUBLE_COMPLEX;
    } else if constexpr (std::is_same_v<T, Boundary<TLattice::NDIM>>) {
        mpi_type = mMPIBoundary;
    }

    assert(mpi_type != MPI_DATATYPE_NULL);

    return mpi_type;
}
#endif
