#pragma once
#include <utility>

#include "../BoundaryModels/BounceBack.hh"
#include "../Collide.hh"
#include "../Data.hh"
#include "../Forces/ChemicalForce.hh"
#include "../Parameters.hh"
#include "ModelBase.hh"

// FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
// Model is given a "TTraits" class that contains stencil, data, force and boundary information

template <class TLattice>
using DefaultTraitFlowField = typename DefaultTrait<TLattice>::template SetBoundary<BounceBack>;

template <class TLattice, class TTraits = DefaultTraitFlowField<TLattice>>
class FlowField : public CollisionBase<TLattice, typename TTraits::Stencil>,
                  public ModelBase<TLattice, TTraits> {  // Inherit from base class to avoid repetition of common
    // calculations

    using Stencil = typename TTraits::Stencil;
    static constexpr int mNDIM = TLattice::NDIM;

   public:
    inline void setTau(double val) {
        mTau = val;
        mInverseTau = 1.0 / mTau;
    }

    inline void collide() override;  // Collision step

    inline void initialise() override;  // Initialisation step

    inline void computeMomenta() override;  // Momenta (density, velocity) calculation

    inline double computeEquilibrium(int k, int idx) override;  // Calculate equilibrium in direction idx

    inline double getTau() { return mTau; }  // Return relaxation time

    template <class, class>
    friend class FlowFieldBinary;

   private:
    double mTau = 1.0;                // TEMPORARY relaxation time
    double mInverseTau = 1.0 / mTau;  // TEMPORARY inverse relaxation time

    std::vector<double>& density = Density<>::get<TLattice>();           // Reference to vector of TDensities
    std::vector<double>& velocity = Velocity<>::get<TLattice, mNDIM>();  // Reference to vector of velocities

    enum { x = 0, y = 1, z = 2 };  // Indices corresponding to x, y, z directions
};

template <class TLattice, class TTraits>
inline void FlowField<TLattice, TTraits>::collide() {  // Collision step

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

        if (this->isCollisionNode(k)) {
            double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

            double equilibriums[Stencil::Q];

            for (int idx = 0; idx < Stencil::Q; idx++) {
                equilibriums[idx] = computeEquilibrium(k, idx);
            }

            this->collisionQ(equilibriums, old_distribution, mInverseTau,
                             k);  // CHANGE NEEDED If no forces, don't require them to be passed
        }
    }
}

template <class TLattice, class TTraits>
inline void FlowField<TLattice, TTraits>::initialise() {  // Initialise model
    this->initialiseProcessors();

    this->mData.generateNeighbors();  // Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this->mt_Forces, mTau, mTau);

#pragma omp parallel for schedule(guided)
    for (int k = 0; k < TLattice::N; k++) {  // loop over k

        Density<>::initialise<TLattice>(1.0, k);  // Set density to 1 initially (This will change)
        Velocity<>::initialise<TLattice, mNDIM>(0.0, k, x);
        if constexpr (mNDIM >= 2) Velocity<>::initialise<TLattice, mNDIM>(0.0, k, y);
        if constexpr (mNDIM == 3) Velocity<>::initialise<TLattice, mNDIM>(0.0, k, z);
    }

    ModelBase<TLattice, TTraits>::mData.communicate(Density<>::getInstance<TLattice>());
    ModelBase<TLattice, TTraits>::mData.communicate(Velocity<>::getInstance<TLattice, TTraits::Lattice::NDIM>());

#pragma omp parallel for schedule(guided)
    for (int k = 0; k < TLattice::N; k++) {  // loop over k

        double* distribution = this->mDistribution.getDistributionPointer(k);
        double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

        for (int idx = 0; idx < Stencil::Q; idx++) {
            double equilibrium = computeEquilibrium(k, idx);
            distribution[idx] = equilibrium;  // Set distributions to equillibrium
            old_distribution[idx] = equilibrium;
        }
    }
}

template <class TLattice, class TTraits>
inline void FlowField<TLattice, TTraits>::computeMomenta() {  // Calculate Density<> and Velocity

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // Loop over k

        if (this->isCollisionNode(k)) {
            double* distribution = this->mDistribution.getDistributionPointer(k);
            velocity[k * Stencil::D + x] =
                this->computeVelocity(distribution, this->mt_Forces, density[k], x, k);  // Calculate velocities

            if constexpr (mNDIM >= 2)
                velocity[k * Stencil::D + y] = this->computeVelocity(distribution, this->mt_Forces, density[k], y, k);
            if constexpr (mNDIM == 3)
                velocity[k * Stencil::D + z] = this->computeVelocity(distribution, this->mt_Forces, density[k], z, k);
            density[k] = this->computeDensity(distribution, k);  // Calculate density
        }
    }
}

template <class TLattice, class TTraits>
inline double FlowField<TLattice, TTraits>::computeEquilibrium(int k, int idx) {
    double gamma = CollisionBase<TLattice, Stencil>::computeGamma(&velocity[k * mNDIM], idx);
    return density[k] * gamma;
}

template <class TMethod>
class PressureForce : public ChemicalForceBinary<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(const int xyz, const int k) {
        return TTraits::Stencil::Cs2 *
                   GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) +
               ChemicalForceBinary<TMethod>::template computeXYZ<TTraits>(xyz, k);
    }
    template <class TTraits>
    inline double computeQ(const int idx, const int k) {
        double sum = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            sum += (TTraits::Stencil::Cs2 *
                        GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) +
                    ChemicalForceBinary<TMethod>::template computeXYZ<TTraits>(xyz, k)) *
                   TTraits::Stencil::Ci_xyz(xyz)[idx];
        }
        return sum;
    }
    template <class TTraits>
    inline double computeDensitySource(int k) {  // SHOULD BE CENTRAL GRADIENTS
        double source = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++)
            source += TTraits::Lattice::DT * 0.5 * TTraits::Stencil::Cs2 *
                      GradientDensity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) *
                      Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz);
        return source;
    }

   private:
};

template <class TLattice, int TNumberOfComponents = 2>
using DefaultTraitFlowFieldPressure =
    typename DefaultTrait<TLattice, TNumberOfComponents>::template SetBoundary<BounceBack>::template SetProcessor<
        Gradients<Density<>, CentralXYZBounceBack>>::template AddForce<PressureForce<He>>;

template <class TLattice, class TTraits = DefaultTraitFlowFieldPressure<TLattice>>
class FlowFieldPressure
    : public CollisionBase<TLattice, typename TTraits::Stencil>,
      public ModelBase<TLattice, TTraits> {  // Inherit from base class to avoid repetition of common
    // calculations

    using Stencil = typename TTraits::Stencil;
    static constexpr int mNDIM = TLattice::NDIM;

   public:
    inline void setTauMin(double val) { mTauMin = val; }
    inline void setTauMax(double val) { mTauMax = val; }

    inline void collide() override;  // Collision step

    inline void initialise() override;  // Initialisation step

    virtual inline void computeMomenta() override;  // Momenta (density, velocity) calculation

    inline double computeEquilibrium(int k, int idx) override;  // Calculate equilibrium in direction idx

    std::vector<double>& pressure = Pressure<>::get<TLattice>();                  // Reference to vector of TDensities
    std::vector<double>& density = Density<>::get<TLattice>();                    // Reference to vector of TDensities
    std::vector<double>& velocity = Velocity<>::get<TLattice, TLattice::NDIM>();  // Reference to vector of velocities

    enum { x = 0, y = 1, z = 2 };  // Indices corresponding to x, y, z directions

    template <class, int, class>
    friend class FlowFieldPressureNComp;

   private:
    double mTauMin = 1;
    double mTauMax = 1;
};

template <class TLattice, class TTraits>
inline void FlowFieldPressure<TLattice, TTraits>::collide() {  // Collision step

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

        if (this->isCollisionNode(k)) {
            double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

            double equilibriums[Stencil::Q];

            for (int idx = 0; idx < Stencil::Q; idx++) {
                equilibriums[idx] = computeEquilibrium(k, idx);
            }

            this->collisionQ(equilibriums, old_distribution, InverseTau<>::get<TLattice>(k),
                             k);  // CHANGE NEEDED If no forces, don't require them to be passed
        }
    }

    ModelBase<TLattice, TTraits>::mData.communicateDistribution();
}

template <class TLattice, class TTraits>
inline void FlowFieldPressure<TLattice, TTraits>::initialise() {  // Initialise model
    this->initialiseProcessors();

    ModelBase<TLattice, TTraits>::mData.generateNeighbors();  // Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this->mt_Forces, mTauMin, mTauMax);

#pragma omp parallel for schedule(guided)
    for (int k = 0; k < TLattice::N; k++) {  // loop over k

        double* distribution = ModelBase<TLattice, TTraits>::mDistribution.getDistributionPointer(k);
        double* old_distribution = ModelBase<TLattice, TTraits>::mDistribution.getDistributionOldPointer(k);
        InverseTau<>::initialise<TLattice>(1.0, k);
        Pressure<>::initialise<TLattice>(1.0, k);  // Set density to 1 initially (This will change)
        Density<>::initialise<TLattice>(1, k);
        Velocity<>::initialise<TLattice, TLattice::NDIM>(0.0, k, x);
        Velocity<>::initialise<TLattice, TLattice::NDIM>(0.0, k, y);
        if constexpr (TLattice::NDIM == 3) Velocity<>::initialise<TLattice, TLattice::NDIM>(0.0, k, z);

        for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {
            double equilibrium = computeEquilibrium(k, idx);
            distribution[idx] = equilibrium;  // Set distributions to equillibrium
            old_distribution[idx] = equilibrium;
        }
    }

    ModelBase<TLattice, TTraits>::mData.communicate(Pressure<>::getInstance<TLattice>());
    ModelBase<TLattice, TTraits>::mData.communicate(Velocity<>::getInstance<TLattice, TLattice::NDIM>());
}

template <class TLattice, class TTraits>
inline void FlowFieldPressure<TLattice, TTraits>::computeMomenta() {  // Calculate Density<> and Velocity

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // Loop over k

        if (this->isCollisionNode(k)) {
            double* distribution = ModelBase<TLattice, TTraits>::mDistribution.getDistributionPointer(k);

            pressure[k] = this->computeDensity(distribution, k);  // Calculate density

            velocity[k * TTraits::Stencil::D + x] =
                1. / (TTraits::Stencil::Cs2) *
                this->computeVelocity(distribution, this->mt_Forces, density[k], x, k);  // Calculate velocities
            velocity[k * TTraits::Stencil::D + y] =
                1. / (TTraits::Stencil::Cs2) * this->computeVelocity(distribution, this->mt_Forces, density[k], y, k);
            if constexpr (TLattice::NDIM == 3)
                velocity[k * TTraits::Stencil::D + z] =
                    1. / (TTraits::Stencil::Cs2) *
                    this->computeVelocity(distribution, this->mt_Forces, density[k], z, k);
        }
    }

    ModelBase<TLattice, TTraits>::mData.communicate(Pressure<>::getInstance<TLattice>());
    ModelBase<TLattice, TTraits>::mData.communicate(Velocity<>::getInstance<TLattice, TLattice::NDIM>());
}

template <class TLattice, class TTraits>
inline double FlowFieldPressure<TLattice, TTraits>::computeEquilibrium(int k, int idx) {
    double velocityFactor = CollisionBase<TLattice, Stencil>::computeVelocityFactor(&velocity[k * mNDIM], idx);
    return Stencil::Weights[idx] * (pressure[k] + density[k] * Stencil::Cs2 * velocityFactor);
}

/** Define a HLBM class for single phase LB simulation that is dervied from the FlowField class.
 * Based on 10.1103/PhysRevE.66.036304, the equation for the equilibrium distribution function needs to be overriden
 * (eq. 10).
 */

template <class TLattice>
using DefaultTraitPorousFlowField = DefaultTraitFlowField<TLattice>;

template <class TLattice, class TTraits = DefaultTraitFlowField<TLattice>>
class PorousFlowField : public FlowField<TLattice, TTraits> {
    using Stencil = typename TTraits::Stencil;
    static constexpr int mNDIM = TLattice::NDIM;

   public:
    inline double computeEquilibrium(int k, int idx) override;

   private:
    static constexpr double mTau = 1.0;                // TEMPORARY relaxation time
    static constexpr double mInverseTau = 1.0 / mTau;  // TEMPORARY inverse relaxation time

    std::vector<double>& density = Density<>::get<TLattice>();           // Reference to vector of TDensities
    std::vector<double>& velocity = Velocity<>::get<TLattice, mNDIM>();  // Reference to vector of velocities
    std::vector<double>& porosity = Porosity<>::get<TLattice>();         // Reference to vector of TPorosities
};

template <class TLattice, class TTraits>
inline double PorousFlowField<TLattice, TTraits>::computeEquilibrium(int k, int idx) {
    // See Eq. 10 in 10.1103/PhysRevE.66.036304

    double velocityFactorFirstOrder =
        CollisionBase<TLattice, Stencil>::computeVelocityFactorFirstOrder(&velocity[k * mNDIM], idx);
    double velocityFactor = CollisionBase<TLattice, Stencil>::computeVelocityFactor(&velocity[k * mNDIM], idx);
    return Stencil::Weights[idx] * density[k] *
           (1 + velocityFactorFirstOrder + 1 / porosity[k] * (velocityFactor - velocityFactorFirstOrder));
}

/* Define a HLBM class for multiphase LB simulation that is dervied from the FlowField class. The main idea
 * is to add a saturation field to the model (instead of Density). The evolution of the saturation field is
 * controlled by the capillary forces.
 */

template <class TLattice>
using DefaultTraitMultiphasePorousFlowField = typename DefaultTrait<TLattice, 2>::template SetBoundary<BounceBack>;

template <class TLattice, class TTraits = DefaultTraitMultiphasePorousFlowField<TLattice>>
class MultiphasePorousFlowField : public FlowField<TLattice, TTraits> {
    using Stencil = typename TTraits::Stencil;
    static constexpr int mNDIM = TLattice::NDIM;

   public:
    inline void initialise() override;  // Initialisation step

    inline void computeMomenta() override;  // Momenta (saturation, velocity) calculation

    inline double computeEquilibrium(int k, int idx) override;  // Calculate equilibrium in direction idx

    inline void setSwIr(double val) { Sw_ir = val; }

   private:
    double mTau = 1.0;                // TEMPORARY relaxation time
    double mInverseTau = 1.0 / mTau;  // TEMPORARY inverse relaxation time

    double Sw_ir = 0.0;

    // Saturation field must be initialised in the main.cc and can be non-unity.
    std::vector<double>& saturation = Saturation<>::get<TLattice>();     // Reference to vector of Saturations
    std::vector<double>& velocity = Velocity<>::get<TLattice, mNDIM>();  // Reference to vector of Velocities
    std::vector<double>& porosity = Porosity<>::get<TLattice>();         // Reference to vector of Porosities

    enum { x = 0, y = 1, z = 2 };  // Indices corresponding to x, y, z directions
};

template <class TLattice, class TTraits>
inline void MultiphasePorousFlowField<TLattice, TTraits>::initialise() {  // Initialise model
    this->initialiseProcessors();

    this->mData.generateNeighbors();  // Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this->mt_Forces, mTau, mTau);

#pragma omp parallel for schedule(guided)
    for (int k = 0; k < TLattice::N; k++) {  // loop over k

        Velocity<>::initialise<TLattice, mNDIM>(0.0, k, x);
        if constexpr (mNDIM >= 2) Velocity<>::initialise<TLattice, mNDIM>(0.0, k, y);
        if constexpr (mNDIM == 3) Velocity<>::initialise<TLattice, mNDIM>(0.0, k, z);
    }

    ModelBase<TLattice, TTraits>::mData.communicate(Velocity<>::getInstance<TLattice, TTraits::Lattice::NDIM>());

#pragma omp parallel for schedule(guided)
    for (int k = 0; k < TLattice::N; k++) {  // loop over k

        double* distribution = this->mDistribution.getDistributionPointer(k);
        double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

        for (int idx = 0; idx < Stencil::Q; idx++) {
            double equilibrium = computeEquilibrium(k, idx);
            distribution[idx] = equilibrium;  // Set distributions to equillibrium
            old_distribution[idx] = equilibrium;
        }
    }
}

template <class TLattice, class TTraits>
inline double MultiphasePorousFlowField<TLattice, TTraits>::computeEquilibrium(int k, int idx) {
    // See Eq. 10 in 10.1103/PhysRevE.66.036304

    double velocityFactorFirstOrder =
        CollisionBase<TLattice, Stencil>::computeVelocityFactorFirstOrder(&velocity[k * mNDIM], idx);
    double velocityFactor = CollisionBase<TLattice, Stencil>::computeVelocityFactor(&velocity[k * mNDIM], idx);
    return Stencil::Weights[idx] * saturation[k] *
           (1 + velocityFactorFirstOrder + 1 / (porosity[k]) * (velocityFactor - velocityFactorFirstOrder));
}

template <class TLattice, class TTraits>
inline void MultiphasePorousFlowField<TLattice, TTraits>::computeMomenta() {  // Calculate  saturation and Velocity

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // Loop over k

        if (this->isCollisionNode(k)) {
            double* distribution = this->mDistribution.getDistributionPointer(k);

            velocity[k * Stencil::D + x] = 0.0;
            if constexpr (mNDIM >= 2) velocity[k * Stencil::D + y] = 0.0;
            if constexpr (mNDIM == 3) velocity[k * Stencil::D + z] = 0.0;

            // if (saturation[k] >= Sw_ir) {
            velocity[k * Stencil::D + x] =
                this->computeVelocity(distribution, this->mt_Forces, saturation[k], x, k);  // Calculate velocities
            if constexpr (mNDIM >= 2)
                velocity[k * Stencil::D + y] =
                    this->computeVelocity(distribution, this->mt_Forces, saturation[k], y, k);
            if constexpr (mNDIM == 3)
                velocity[k * Stencil::D + z] =
                    this->computeVelocity(distribution, this->mt_Forces, saturation[k], z, k);
            // }
            saturation[k] = this->computeDensity(distribution, k);  // Calculate saturation
        }
    }
}

// The main difference between the MultiphasePorousFlowField and MultiphasePorousFlowFieldInSalt
// is that the former uses the BounceBack boundary condition, while the latter uses the BounceBackHLBM
// boundary condition. Also, for the velocity calculation, the former uses Velocity<>::get<TLattice, mNDIM>()
// while the latter uses VelocityPorous<>::get<TLattice, mNDIM>().

template <class TLattice>
using DefaultTraitMultiphasePorousFlowFieldInSalt = typename DefaultTrait<TLattice, 2>::template SetBoundary<BounceBack>;

template <class TLattice, class TTraits = DefaultTraitMultiphasePorousFlowField<TLattice>>
class MultiphasePorousFlowFieldInSalt : public FlowField<TLattice, TTraits> {
    using Stencil = typename TTraits::Stencil;
    static constexpr int mNDIM = TLattice::NDIM;

   public:
    inline void initialise() override;  // Initialisation step

    inline void computeMomenta() override;  // Momenta (saturation, velocityPorous) calculation

    inline double computeEquilibrium(int k, int idx) override;  // Calculate equilibrium in direction idx

    inline void setSwIr(double val) { Sw_ir = val; }

   private:
    double mTau = 1.0;                // TEMPORARY relaxation time
    double mInverseTau = 1.0 / mTau;  // TEMPORARY inverse relaxation time

    double Sw_ir = 0.0;

    // Saturation field must be initialised in the main.cc and can be non-unity.
    std::vector<double>& saturation = Saturation<>::get<TLattice>();  // Reference to vector of Saturations
    std::vector<double>& velocityPorous =
        VelocityPorous<>::get<TLattice, mNDIM>();                 // Reference to vector of Velocities
    std::vector<double>& porosity = Porosity<>::get<TLattice>();  // Reference to vector of Porosities

    enum { x = 0, y = 1, z = 2 };  // Indices corresponding to x, y, z directions
};

template <class TLattice, class TTraits>
inline void MultiphasePorousFlowFieldInSalt<TLattice, TTraits>::initialise() {  // Initialise model
    this->initialiseProcessors();

    this->mData.generateNeighbors();  // Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this->mt_Forces, mTau, mTau);

#pragma omp parallel for schedule(guided)
    for (int k = 0; k < TLattice::N; k++) {  // loop over k

        VelocityPorous<>::initialise<TLattice, mNDIM>(0.0, k, x);
        if constexpr (mNDIM >= 2) VelocityPorous<>::initialise<TLattice, mNDIM>(0.0, k, y);
        if constexpr (mNDIM == 3) VelocityPorous<>::initialise<TLattice, mNDIM>(0.0, k, z);
    }

    ModelBase<TLattice, TTraits>::mData.communicate(VelocityPorous<>::getInstance<TLattice, TTraits::Lattice::NDIM>());

#pragma omp parallel for schedule(guided)
    for (int k = 0; k < TLattice::N; k++) {  // loop over k

        double* distribution = this->mDistribution.getDistributionPointer(k);
        double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

        for (int idx = 0; idx < Stencil::Q; idx++) {
            double equilibrium = computeEquilibrium(k, idx);
            distribution[idx] = equilibrium;  // Set distributions to equillibrium
            old_distribution[idx] = equilibrium;
        }
    }
}

template <class TLattice, class TTraits>
inline double MultiphasePorousFlowFieldInSalt<TLattice, TTraits>::computeEquilibrium(int k, int idx) {
    // See Eq. 10 in 10.1103/PhysRevE.66.036304

    double velocityFactorFirstOrder =
        CollisionBase<TLattice, Stencil>::computeVelocityFactorFirstOrder(&velocityPorous[k * mNDIM], idx);
    double velocityFactor = CollisionBase<TLattice, Stencil>::computeVelocityFactor(&velocityPorous[k * mNDIM], idx);
    return Stencil::Weights[idx] * saturation[k] *
           (1 + velocityFactorFirstOrder + 1 / (porosity[k]) * (velocityFactor - velocityFactorFirstOrder));
}

template <class TLattice, class TTraits>
inline void
MultiphasePorousFlowFieldInSalt<TLattice, TTraits>::computeMomenta() {  // Calculate  saturation and velocityPorous

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // Loop over k

        if (this->isCollisionNode(k)) {
            double* distribution = this->mDistribution.getDistributionPointer(k);

            velocityPorous[k * Stencil::D + x] = 0.0;
            if constexpr (mNDIM >= 2) velocityPorous[k * Stencil::D + y] = 0.0;
            if constexpr (mNDIM == 3) velocityPorous[k * Stencil::D + z] = 0.0;

            // if (saturation[k] >= Sw_ir) {
            velocityPorous[k * Stencil::D + x] =
                this->computeVelocity(distribution, this->mt_Forces, saturation[k], x, k);  // Calculate velocities
            if constexpr (mNDIM >= 2)
                velocityPorous[k * Stencil::D + y] =
                    this->computeVelocity(distribution, this->mt_Forces, saturation[k], y, k);
            if constexpr (mNDIM == 3)
                velocityPorous[k * Stencil::D + z] =
                    this->computeVelocity(distribution, this->mt_Forces, saturation[k], z, k);
            // }
            saturation[k] = this->computeDensity(distribution, k);  // Calculate saturation
        }
    }
}