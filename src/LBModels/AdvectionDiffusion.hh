#pragma once
#include <utility>

#include "../AddOns/AddOns.hh"
#include "../BoundaryModels/Boundaries.hh"
#include "../Collide.hh"
#include "../Data.hh"
#include "../Forces/Forces.hh"
#include "../GradientStencils/GradientStencils.hh"
#include "../Parallel.hh"
#include "../Parameters.hh"
#include "ModelBase.hh"

// FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
// Model is given a "TTraits" class that contains stencil, data, force and boundary information

// Needs to be adapted for multicomponent zeroth moment vals

template <class TLattice>
using DefaultTraitAdvectionDiffusion = DefaultTrait<TLattice>;

template <class TZerothMoment, class TLattice, class TTraits = DefaultTraitAdvectionDiffusion<TLattice>>
class AdvectionDiffusion
    : public CollisionBase<TLattice, typename TTraits::Stencil>,
      public ModelBase<TLattice, TTraits> {  // Inherit from base class to avoid repetition of common
    // calculations

    using Stencil = typename TTraits::Stencil;
    static constexpr int mNDIM = TLattice::NDIM;

   public:
    inline void collide() override;  // Collision step

    inline void initialise() override;  // Initialisation step

    inline void computeMomenta() override;  // Momenta (density, velocity) calculation

    inline double computeEquilibrium(int k, int idx) override;  // Calculate equilibrium in direction idx

    inline void setDiffusivity(double D) {
        mTau = D / Stencil::Cs2 + 0.5 * TLattice::DT;
        mInverseTau = 1.0 / mTau;
    }

   private:
    double mTau = 0.02 / Stencil::Cs2 + 0.5 * TLattice::DT;  // TEMPORARY relaxation time
    double mInverseTau = 1.0 / mTau;                         // TEMPORARY inverse relaxation time

    std::vector<double>& zerothmoment = TZerothMoment::template get<TLattice>();  // Reference to vector of TDensities
    std::vector<double>& velocity = Velocity<>::get<TLattice, mNDIM>();           // Reference to vector of velocities

    enum { x = 0, y = 1, z = 2 };  // Indices corresponding to x, y, z directions
};

template <class TZerothMoment, class TLattice, class TTraits>
inline void AdvectionDiffusion<TZerothMoment, TLattice, TTraits>::collide() {  // Collision step

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

    this->mData.communicateDistribution();
}

template <class TZerothMoment, class TLattice, class TTraits>
inline void AdvectionDiffusion<TZerothMoment, TLattice, TTraits>::initialise() {  // Initialise model
    this->initialiseProcessors();

    this->mData.generateNeighbors();  // Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this->mt_Forces, mTau, mTau);

#pragma omp parallel for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

        TZerothMoment::template initialise<TLattice>(0.0, k);  // Set density to 1 initially (This will change)
    }

    this->mData.communicate(BoundaryLabels<TLattice::NDIM>::template getInstance<TLattice>());
    this->mData.communicate(TZerothMoment::template getInstance<TLattice>());

#pragma omp parallel for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {
        double* distribution = this->mDistribution.getDistributionPointer(k);
        double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

        for (int idx = 0; idx < Stencil::Q; idx++) {
            double equilibrium = computeEquilibrium(k, idx);
            distribution[idx] = equilibrium;  // Set distributions to equillibrium
            old_distribution[idx] = equilibrium;
        }
    }
}

template <class TZerothMoment, class TLattice, class TTraits>
inline void
AdvectionDiffusion<TZerothMoment, TLattice, TTraits>::computeMomenta() {  // Calculate Density<> and Velocity

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // Loop over k

        if (this->isCollisionNode(k)) {
            double* distribution = this->mDistribution.getDistributionPointer(k);

            zerothmoment[k] = this->computeDensity(distribution, k);  // Calculate density
        }
    }
}

template <class TZerothMoment, class TLattice, class TTraits>
inline double AdvectionDiffusion<TZerothMoment, TLattice, TTraits>::computeEquilibrium(int k, int idx) {
    double gamma = CollisionBase<TLattice, Stencil>::computeGamma(&velocity[k * mNDIM], idx);
    return zerothmoment[k] * gamma;
}

template <class TZerothMoment, class TZerothMomentOld, class TMethod = AdvectionDiffusionSourceMethod>
class AdvectionDiffusionCorrection : public ForceBase<TMethod> {
   public:
    template <class TTraits>
    inline double computeXYZ(int xyz, int k) const {
        return (TZerothMoment::template get<typename TTraits::Lattice>(k) *
                    Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz) -
                TZerothMomentOld::template get<typename TTraits::Lattice>(k) *
                    VelocityOld<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, xyz)) /
               TTraits::Lattice::DT;
    }
};

template <class TLattice, class TZerothMoment, class TZerothMomentOld>
using DefaultTraitAdvectionDiffusionCorrected =
    typename DefaultTrait<TLattice>::template SetForce<AdvectionDiffusionCorrection<TZerothMoment, TZerothMomentOld>>;

template <class TZerothMoment, class TZerothMomentOld, class TLattice,
          class TTraits = DefaultTraitAdvectionDiffusionCorrected<TLattice, TZerothMoment, TZerothMomentOld>>
class AdvectionDiffusionCorrected
    : public CollisionBase<TLattice, typename TTraits::Stencil>,
      public ModelBase<TLattice, TTraits> {  // Inherit from base class to avoid repetition of common
    // calculations

    using Stencil = typename TTraits::Stencil;
    static constexpr int mNDIM = TLattice::NDIM;

   public:
    inline void collide() override;  // Collision step

    inline void initialise() override;  // Initialisation step

    inline void computeMomenta() override;  // Momenta (density, velocity) calculation

    inline double computeEquilibrium(int k, int idx) override;  // Calculate equilibrium in direction idx

    inline void setDiffusivity(double D) {
        mTau = D / Stencil::Cs2 + 0.5 * TLattice::DT;
        mInverseTau = 1.0 / mTau;
    }

   private:
    double mTau = 0.02 / Stencil::Cs2 + 0.5 * TLattice::DT;  // TEMPORARY relaxation time
    double mInverseTau = 1.0 / mTau;                         // TEMPORARY inverse relaxation time

    std::vector<double>& zerothmoment = TZerothMoment::template get<TLattice>();  // Reference to vector of TDensities
    std::vector<double>& zerothmomentold =
        TZerothMomentOld::template get<TLattice>();                      // Reference to vector of TDensities
    std::vector<double>& velocity = Velocity<>::get<TLattice, mNDIM>();  // Reference to vector of velocities

    enum { x = 0, y = 1, z = 2 };  // Indices corresponding to x, y, z directions
};

template <class TZerothMoment, class TZerothMomentOld, class TLattice, class TTraits>
inline void
AdvectionDiffusionCorrected<TZerothMoment, TZerothMomentOld, TLattice, TTraits>::collide() {  // Collision step

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

    this->mData.communicateDistribution();
}

template <class TZerothMoment, class TZerothMomentOld, class TLattice, class TTraits>
inline void
AdvectionDiffusionCorrected<TZerothMoment, TZerothMomentOld, TLattice, TTraits>::initialise() {  // Initialise model
    this->initialiseProcessors();

    this->mData.generateNeighbors();  // Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this->mt_Forces, mTau, mTau);

#pragma omp parallel for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

        TZerothMoment::template initialise<TLattice>(0.0, k);  // Set density to 1 initially (This will change)
        TZerothMomentOld::template initialise<TLattice>(TZerothMoment::template get<TLattice>(k),
                                                        k);  // Set density to 1 initially (This will change)
    }

    this->mData.communicate(BoundaryLabels<TLattice::NDIM>::template getInstance<TLattice>());
    this->mData.communicate(TZerothMoment::template getInstance<TLattice>());
    this->mData.communicate(TZerothMomentOld::template getInstance<TLattice>());

#pragma omp parallel for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {
        double* distribution = this->mDistribution.getDistributionPointer(k);
        double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

        for (int idx = 0; idx < Stencil::Q; idx++) {
            double equilibrium = computeEquilibrium(k, idx);
            distribution[idx] = equilibrium;  // Set distributions to equillibrium
            old_distribution[idx] = equilibrium;
        }

        // std::cout<<Geometry<TLattice>::getBoundaryType(k)<<std::endl;
    }
}

template <class TZerothMoment, class TZerothMomentOld, class TLattice, class TTraits>
inline void AdvectionDiffusionCorrected<TZerothMoment, TZerothMomentOld, TLattice,
                                        TTraits>::computeMomenta() {  // Calculate Density<> and Velocity

    TZerothMoment::template get<typename TTraits::Lattice>().swap(
        TZerothMomentOld::template get<typename TTraits::Lattice>());

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // Loop over k

        if (this->isCollisionNode(k)) {
            double* distribution = this->mDistribution.getDistributionPointer(k);

            zerothmoment[k] = this->computeDensity(distribution, k);  // Calculate density
            // humidity[k] = this -> computeDensity(distribution, k); //Calculate density
        }
    }

    this->mData.communicate(TZerothMomentOld::template getInstance<TLattice>());
    this->mData.communicate(TZerothMoment::template getInstance<TLattice>());
}

template <class TZerothMoment, class TZerothMomentOld, class TLattice, class TTraits>
inline double AdvectionDiffusionCorrected<TZerothMoment, TZerothMomentOld, TLattice, TTraits>::computeEquilibrium(
    int k, int idx) {
    double gamma = CollisionBase<TLattice, Stencil>::computeGamma(&velocity[k * mNDIM], idx);
    return zerothmoment[k] * gamma;
}
