#pragma once
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Collide.hh"
#include "Geometry.hh"
#include "Global.hh"
#include "Lattice.hh"
#include "Service.hh"
#include "Stencil.hh"

/**
 * \file Forcing.hh
 * \brief Contains commonly used forcing/source models for the LBM.
 * There is a base class which contains the dimensions of each force needed and the dependence on the relaxation time.
 */

template <class... TStencils>
struct ForcingBase {
    using mt_Stencils = std::tuple<TStencils...>;
    using Prefactor = NoTauDependence;

    inline void reset() {}

    template <class TTraits, class TStencil>
    inline static constexpr int StencilToInt() {
        if constexpr (std::is_same_v<TStencil, Cartesian>)
            return TTraits::Lattice::NDIM;
        else if constexpr (std::is_same_v<TStencil, AllDirections>)
            return TTraits::Stencil::Q;
        else if constexpr (std::is_same_v<TStencil, One>)
            return 1;
        else
            return TTraits::Lattice::NDIM;
    }

    template <class TTraits>
    inline static constexpr auto ForceStencils() -> const std::array<size_t, sizeof...(TStencils)> {
        constexpr std::array<size_t, sizeof...(TStencils)> Stencilarray = {StencilToInt<TTraits, TStencils>()...};
        return Stencilarray;
    }
};

struct SimpleForcing : ForcingBase<Cartesian> {
    std::vector<double> ma_Force;

    using Prefactor = GuoPrefactor;

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        if (ma_Force.size() < TTraits::Lattice::NDIM) ma_Force.resize(TTraits::Lattice::NDIM, 0);
        ma_Force[0] += f.template computeXYZ<TTraits>(0, k);
        if constexpr (TTraits::Lattice::NDIM >= 2) ma_Force[1] += f.template computeXYZ<TTraits>(1, k);
        if constexpr (TTraits::Lattice::NDIM == 3) ma_Force[2] += f.template computeXYZ<TTraits>(2, k);
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing

        double prefactor = TTraits::Stencil::Weights[idx];  // Prefactor for Guo forcing

        double ci_dot_force = (TTraits::Stencil::Ci_x[idx] * ma_Force[0]);
        if constexpr (TTraits::Stencil::D > 1) ci_dot_force += (TTraits::Stencil::Ci_y[idx] * ma_Force[1]);
        if constexpr (TTraits::Stencil::D > 2) ci_dot_force += (TTraits::Stencil::Ci_z[idx] * ma_Force[2]);

        return prefactor * ci_dot_force / TTraits::Stencil::Cs2;
    }
};

struct ADSinkSourceMethod : ForcingBase<Cartesian> {
    double m_Source0D;

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        m_Source0D += f.template compute<TTraits>(k);
    }

    template <class TTraits>
    inline double compute(int idx, int k) {
        double prefactor = TTraits::Stencil::Weights[idx];

        return prefactor * m_Source0D;
    }
};

struct Guo : ForcingBase<Cartesian> {
    std::vector<double> ma_Force;

    using Prefactor = GuoPrefactor;

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        if (ma_Force.size() < TTraits::Lattice::NDIM) ma_Force.resize(TTraits::Lattice::NDIM, 0);
        ma_Force[0] += f.template computeXYZ<TTraits>(0, k);
        if constexpr (TTraits::Lattice::NDIM >= 2) ma_Force[1] += f.template computeXYZ<TTraits>(1, k);
        if constexpr (TTraits::Lattice::NDIM == 3) ma_Force[2] += f.template computeXYZ<TTraits>(2, k);
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing

        double prefactor = TTraits::Stencil::Weights[idx];  // Prefactor for Guo forcing

        double ci_dot_velocity =
            (TTraits::Stencil::Ci_x[idx] * Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0));
        if constexpr (TTraits::Stencil::D > 1)
            ci_dot_velocity += (TTraits::Stencil::Ci_y[idx] *
                                Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 1));
        if constexpr (TTraits::Stencil::D > 2)
            ci_dot_velocity += (TTraits::Stencil::Ci_z[idx] *
                                Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 2));

        double forceterm =
            prefactor *
            (((TTraits::Stencil::Ci_x[idx] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0)) /
                  TTraits::Stencil::Cs2 +
              ci_dot_velocity * TTraits::Stencil::Ci_x[idx] / (TTraits::Stencil::Cs2 * TTraits::Stencil::Cs2)) *
             ma_Force[0]);  // Force Calculation
        if constexpr (TTraits::Stencil::D > 1)
            forceterm +=
                prefactor *
                (((TTraits::Stencil::Ci_y[idx] -
                   Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 1)) /
                      TTraits::Stencil::Cs2 +
                  ci_dot_velocity * TTraits::Stencil::Ci_y[idx] / (TTraits::Stencil::Cs2 * TTraits::Stencil::Cs2)) *
                 ma_Force[1]);
        if constexpr (TTraits::Stencil::D > 2)
            forceterm +=
                prefactor *
                (((TTraits::Stencil::Ci_z[idx] -
                   Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 2)) /
                      TTraits::Stencil::Cs2 +
                  ci_dot_velocity * TTraits::Stencil::Ci_z[idx] / (TTraits::Stencil::Cs2 * TTraits::Stencil::Cs2)) *
                 ma_Force[2]);

        return forceterm;
    }
};

struct WellBalancedForce : ForcingBase<Cartesian> {
    Guo mGuo;

    using Prefactor = GuoPrefactor;

    double mForceDotVelocity = 0;

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        mGuo.precompute<TTraits>(f, k);
        mForceDotVelocity = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++)
            mForceDotVelocity += GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) *
                                 Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz);
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing

        double ci_dot_velocity = 0;
        double ci_dot_ci = 0;
        double forceterm = 0;

        double prefactor = TTraits::Stencil::Weights[idx];  // Prefactor for Guo forcing

        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            ci_dot_velocity += (TTraits::Stencil::Ci_xyz(xyz)[idx] *
                                Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(
                                    k, xyz));  // Dot product of discrete velocity vector
                                               // with velocity
            ci_dot_ci += TTraits::Stencil::Ci_xyz(xyz)[idx] * TTraits::Stencil::Ci_xyz(xyz)[idx];
        }

        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            forceterm += GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) *
                         ci_dot_velocity * TTraits::Stencil::Ci_xyz(xyz)[idx] / TTraits::Stencil::Cs2;  // Force
                                                                                                        // Calculation
        }

        return mGuo.compute<TTraits>(idx, k) + prefactor * forceterm +
               prefactor * (-mForceDotVelocity +
                            0.5 * ((ci_dot_ci) / TTraits::Stencil::Cs2 - TTraits::Lattice::NDIM) * mForceDotVelocity);
    }
};

struct AllenCahnSourceMethod : ForcingBase<Cartesian> {
    std::vector<double> ma_Force;

    using Prefactor = GuoPrefactor;

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        if (ma_Force.size() < TTraits::Lattice::NDIM) ma_Force.resize(TTraits::Lattice::NDIM, 0);
        ma_Force[0] += f.template computeXYZ<TTraits>(0, k);
        if constexpr (TTraits::Lattice::NDIM >= 2) ma_Force[1] += f.template computeXYZ<TTraits>(1, k);
        if constexpr (TTraits::Lattice::NDIM == 3) ma_Force[2] += f.template computeXYZ<TTraits>(2, k);
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing

        double forceterm = 0;

        double prefactor = TTraits::Stencil::Weights[idx];  // Prefactor for Guo forcing

        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            forceterm += prefactor * (TTraits::Stencil::Ci_xyz(xyz)[idx]) * ma_Force[xyz];  // Force
                                                                                            // Calculation
        }

        return forceterm / (TTraits::Stencil::Cs2 * TTraits::Lattice::DT * TTraits::Lattice::DT);
    }
};

struct AdvectionDiffusionSourceMethod : ForcingBase<Cartesian> {
    std::vector<double> ma_Source;

    using Prefactor = GuoPrefactor;

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        if (ma_Source.size() < TTraits::Lattice::NDIM) ma_Source.resize(TTraits::Lattice::NDIM, 0);
        ma_Source[0] += f.template computeXYZ<TTraits>(0, k);
        if constexpr (TTraits::Lattice::NDIM >= 2) ma_Source[1] += f.template computeXYZ<TTraits>(1, k);
        if constexpr (TTraits::Lattice::NDIM == 3) ma_Source[2] += f.template computeXYZ<TTraits>(2, k);
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing

        double sourceterm = 0;

        double prefactor = TTraits::Stencil::Weights[idx];  // Prefactor for Guo forcing

        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            sourceterm += (TTraits::Stencil::Ci_xyz(xyz)[idx]) * ma_Source[xyz];  // Force
                                                                                  // Calculation
        }

        return prefactor * (sourceterm / (TTraits::Stencil::Cs2));
    }
};

struct EvaporationSourceMethod : ForcingBase<Cartesian, One> {
    std::vector<double> ma_Source;
    double m_Source0D;

    using Prefactor = GuoPrefactor;

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        m_Source0D += f.template compute<TTraits>(k);
        if (ma_Source.size() < TTraits::Lattice::NDIM) ma_Source.resize(TTraits::Lattice::NDIM, 0);
        ma_Source[0] += f.template computeXYZ<TTraits>(0, k);
        if constexpr (TTraits::Lattice::NDIM >= 2) ma_Source[1] += f.template computeXYZ<TTraits>(1, k);
        if constexpr (TTraits::Lattice::NDIM == 3) ma_Source[2] += f.template computeXYZ<TTraits>(2, k);
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing

        double sourceterm = 0;

        double prefactor = TTraits::Stencil::Weights[idx];  // Prefactor for Guo forcing

        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            sourceterm += (TTraits::Stencil::Ci_xyz(xyz)[idx]) * ma_Source[xyz];  // Force
                                                                                  // Calculation
        }

        return prefactor * (m_Source0D + sourceterm / (TTraits::Stencil::Cs2));
    }
};

struct He : ForcingBase<Cartesian> {
    std::vector<double> ma_Force;

    using Prefactor = GuoPrefactor;

    double velocity_dot_force = 0;

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        if (ma_Force.size() < TTraits::Lattice::NDIM) ma_Force.resize(TTraits::Lattice::NDIM, 0);
        ma_Force[0] += f.template computeXYZ<TTraits>(0, k);
        if constexpr (TTraits::Lattice::NDIM >= 2) ma_Force[1] += f.template computeXYZ<TTraits>(1, k);
        if constexpr (TTraits::Lattice::NDIM == 3) ma_Force[2] += f.template computeXYZ<TTraits>(2, k);

        velocity_dot_force = 0;
        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            velocity_dot_force +=
                (ma_Force[xyz] * Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(
                                     k, xyz));  // Dot product of discrete velocity vector with velocity
        }
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing

        double prefactor = CollisionBase<typename TTraits::Lattice, typename TTraits::Stencil>::computeGamma(
            &Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0),
            idx);  // Prefactor for Guo forcing

        double ci_dot_force = 0;

        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            ci_dot_force += (TTraits::Stencil::Ci_xyz(xyz)[idx] * ma_Force[xyz]);  // Dot product of discrete velocity
                                                                                   // vector with velocity
        }

        return prefactor * (ci_dot_force - velocity_dot_force);
    }
};

struct NCompForce : ForcingBase<Cartesian> {
    He mHe;

    using Prefactor = GuoPrefactor;

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        mHe.precompute<TTraits>(f, k);
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing

        double forceterm = 0;

        double prefactor = TTraits::Stencil::Weights[idx];  // Prefactor for Guo forcing

        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            forceterm += (TTraits::Stencil::Ci_xyz(xyz)[idx] -
                          Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz)) *
                         GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) *
                         TTraits::Stencil::Cs2 *
                         (CollisionBase<typename TTraits::Lattice, typename TTraits::Stencil>::computeGamma(
                              &Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0), idx) -
                          prefactor);  // Force
                                       // Calculation
        }

        return mHe.compute<TTraits>(idx, k) + forceterm;
    }
};

struct Lee : ForcingBase<Cartesian, AllDirections> {
    int kPrecompute;
    std::vector<double> ma_Force;
    std::vector<double> ma_ForceQ;

    double velocity_dot_force = 0;

    inline void reset() {
        std::fill(ma_Force.begin(), ma_Force.end(), 0);
        std::fill(ma_ForceQ.begin(), ma_ForceQ.end(), 0);
    }

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        if (ma_Force.size() < TTraits::Lattice::NDIM) ma_Force.resize(TTraits::Lattice::NDIM, 0);
        if (ma_ForceQ.size() < TTraits::Stencil::Q) ma_ForceQ.resize(TTraits::Stencil::Q, 0);
        ma_Force[0] += f.template computeXYZ<TTraits>(0, k);
        if constexpr (TTraits::Lattice::NDIM >= 2) ma_Force[1] += f.template computeXYZ<TTraits>(1, k);
        if constexpr (TTraits::Lattice::NDIM == 3) ma_Force[2] += f.template computeXYZ<TTraits>(2, k);
        for (int q = 0; q < TTraits::Stencil::Q; q++) ma_ForceQ[q] += f.template computeQ<TTraits>(q, k);

        velocity_dot_force = 0;
        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            velocity_dot_force +=
                (ma_Force[xyz] * Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(
                                     k, xyz));  // Dot product of discrete velocity vector with velocity
        }
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing

        double prefactor =
            TTraits::Lattice::DT * CollisionBase<typename TTraits::Lattice, typename TTraits::Stencil>::computeGamma(
                                       &Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0),
                                       idx);  // Prefactor for Guo forcing

        return prefactor * (ma_ForceQ[idx] - velocity_dot_force);
    }
};

struct LeeEquilibrium : Lee {
    using Prefactor = SplitGuoPrefactor;
};

struct LeeGamma0 : ForcingBase<Cartesian, AllDirections> {
    std::vector<double> ma_Force;
    std::vector<double> ma_ForceQ;

    double velocity_dot_force = 0;

    inline void reset() {
        std::fill(ma_Force.begin(), ma_Force.end(), 0);
        std::fill(ma_ForceQ.begin(), ma_ForceQ.end(), 0);
    }

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        if (ma_Force.size() < TTraits::Lattice::NDIM) ma_Force.resize(TTraits::Lattice::NDIM, 0);
        if (ma_ForceQ.size() < TTraits::Stencil::Q) ma_ForceQ.resize(TTraits::Stencil::Q, 0);
        ma_Force[0] += f.template computeXYZ<TTraits>(0, k);
        if constexpr (TTraits::Lattice::NDIM >= 2) ma_Force[1] += f.template computeXYZ<TTraits>(1, k);
        if constexpr (TTraits::Lattice::NDIM == 3) ma_Force[2] += f.template computeXYZ<TTraits>(2, k);
        for (int q = 0; q < TTraits::Stencil::Q; q++) ma_ForceQ[q] += f.template computeQ<TTraits>(q, k);

        velocity_dot_force = 0;
        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            velocity_dot_force +=
                (ma_Force[xyz] * Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(
                                     k, xyz));  // Dot product of discrete velocity vector with velocity
        }
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing

        double prefactor = TTraits::Lattice::DT * TTraits::Stencil::Weights[idx];  // Prefactor for Guo forcing

        return prefactor * (ma_ForceQ[idx] - velocity_dot_force);
    }
};

struct LeeEquilibriumGamma0 : LeeGamma0 {
    using Prefactor = SplitGuoPrefactor;
};

struct HeGamma0 : ForcingBase<Cartesian> {
    std::vector<double> ma_Force;

    double velocity_dot_force = 0;
    double ci_dot_force = 0;

    inline void reset() { std::fill(ma_Force.begin(), ma_Force.end(), 0); }

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int k) {
        if (ma_Force.size() < TTraits::Lattice::NDIM) ma_Force.resize(TTraits::Lattice::NDIM, 0);
        ma_Force[0] += f.template computeXYZ<TTraits>(0, k);
        if constexpr (TTraits::Lattice::NDIM >= 2) ma_Force[1] += f.template computeXYZ<TTraits>(1, k);
        if constexpr (TTraits::Lattice::NDIM == 3) ma_Force[2] += f.template computeXYZ<TTraits>(2, k);

        velocity_dot_force = 0;
        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            velocity_dot_force +=
                (ma_Force[xyz] * Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(
                                     k, xyz));  // Dot product of discrete velocity vector with velocity
        }
        for (int xyz = 0; xyz < TTraits::Stencil::D; xyz++) {
            ci_dot_force += (ma_Force[xyz] * Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(
                                                 k, xyz));  // Dot product of discrete velocity vector with velocity
        }
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing

        double prefactor = TTraits::Lattice::DT * TTraits::Stencil::Weights[idx];  // Prefactor for Guo forcing

        return prefactor * (ci_dot_force - velocity_dot_force);
    }
};

/*
struct LeeMuLocal : ForcingBase<AllDirections> {

    std::vector<double> ma_ForceQ;

    template<class TTraits, class TForce>
    const inline void precompute(TForce& f, int k){

        if (ma_ForceQ.size()<TTraits::Stencil::Q) ma_ForceQ.resize(TTraits::Stencil::Q,0);
        for (int q = 0; q < TTraits::Stencil::Q; q++) ma_ForceQ[q] += f.template computeQ<TTraits>(q,k);

    }

    template<class TTraits, class TForce>
    inline void precompute(TForce& f, int q, int k){

        if (ma_ForceQ.size()<TTraits::Stencil::Q) ma_ForceQ.resize(TTraits::Stencil::Q,0);
        ma_ForceQ[q] += f.template computeQ<TTraits>(q,k);

    }

    template<class TTraits>
    inline double compute(int idx, int k) { //Guo forcing

        double prefactor = 0.5 * TTraits::Lattice::DT * CollisionBase<typename TTraits::Lattice,typename
TTraits::Stencil>::computeGamma(&Velocity<>::get<typename TTraits::Lattice,TTraits::Lattice::NDIM>(k,0),idx);
//Prefactor for Guo forcing

        return prefactor*(ma_ForceQ[idx]);

    }
};


struct LeeMuNonLocal : ForcingBase<AllDirections> {

    std::vector<double> ma_ForceQ;

    template<class TTraits, class TForce>
    const inline void precompute(TForce& f, int k){

        if (ma_ForceQ.size()<TTraits::Stencil::Q) ma_ForceQ.resize(TTraits::Stencil::Q,0);
        for (int q = 0; q < TTraits::Stencil::Q; q++) ma_ForceQ[q] += f.template computeQ<TTraits>(q,k);

    }

    template<class TTraits, class TForce>
    const inline void precompute(TForce& f, int q, int k){

        if (ma_ForceQ.size()<TTraits::Stencil::Q) ma_ForceQ.resize(TTraits::Stencil::Q,0);
        ma_ForceQ[q] += f.template computeQ<TTraits>(q,k);

    }

    template<class TTraits>
    const inline double compute(int idx, int k) { //Guo forcing

        using data = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;

        double gamma;

        if (Geometry<typename TTraits::Lattice>::getBoundaryType(data::getInstance().getNeighbors()[k *
TTraits::Stencil::Q + idx]) != 0) gamma = CollisionBase<typename TTraits::Lattice,typename
TTraits::Stencil>::computeGamma(&Velocity<>::get<typename
TTraits::Lattice,TTraits::Lattice::NDIM>(k,0),TTraits::Stencil::Opposites[idx]); else gamma = CollisionBase<typename
TTraits::Lattice,typename TTraits::Stencil>::computeGamma(&Velocity<>::get<typename
TTraits::Lattice,TTraits::Lattice::NDIM>(data::getInstance().getNeighbor(k,idx),0),idx);

        const double prefactor = 0.5 * TTraits::Lattice::DT * gamma; //Prefactor for Guo forcing

        return prefactor*(ma_ForceQ[idx]);

    }
};
*/

struct LeeMuLocal : ForcingBase<AllDirections> {
    std::vector<double> ma_ForceQ;
    template <class TTraits, class TForce>
    const inline void precompute(TForce& f, int k) {
        ma_ForceQ = {};
        for (int q = 0; q < TTraits::Stencil::Q; q++) ma_ForceQ.push_back(f.template computeQ<TTraits>(q, k));
    }

    template <class TTraits, class TForce>
    inline void precompute(TForce& f, int q, int k) {
        ma_ForceQ.push_back(f.template computeQ<TTraits>(q, k));
    }

    template <class TTraits>
    inline double compute(int idx, int k) {  // Guo forcing
        double prefactor = 0.5 * TTraits::Lattice::DT *
                           CollisionBase<typename TTraits::Lattice, typename TTraits::Stencil>::computeGamma(
                               &Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0),
                               idx);  // Prefactor for Guo forcing
        return prefactor * (ma_ForceQ[idx]);
    }
};
struct LeeMuNonLocal : ForcingBase<AllDirections> {
    std::vector<double> ma_ForceQ;
    template <class TTraits, class TForce>
    const inline void precompute(TForce& f, int k) {
        ma_ForceQ = {};
        for (int q = 0; q < TTraits::Stencil::Q; q++) ma_ForceQ.push_back(f.template computeQ<TTraits>(q, k));
    }

    template <class TTraits, class TForce>
    const inline void precompute(TForce& f, int q, int k) {
        ma_ForceQ.push_back(f.template computeQ<TTraits>(q, k));
    }

    template <class TTraits>
    const inline double compute(int idx, int k) {  // Guo forcing
        using data = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;
        double gamma;
        if (!Geometry<typename TTraits::Lattice>::isBoundary(
                data::getInstance().getNeighbors()[k * TTraits::Stencil::Q + idx]))
            gamma = CollisionBase<typename TTraits::Lattice, typename TTraits::Stencil>::computeGamma(
                &Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(
                    data::getInstance().getNeighbors()[k * TTraits::Stencil::Q + idx], 0),
                idx);
        else
            gamma = CollisionBase<typename TTraits::Lattice, typename TTraits::Stencil>::computeGamma(
                &Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0), idx);
        const double prefactor = 0.5 * TTraits::Lattice::DT * gamma;  // Prefactor for Guo forcing
        return prefactor * (ma_ForceQ[idx]);
    }
};