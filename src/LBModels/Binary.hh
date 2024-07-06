#pragma once
#include <omp.h>

#include <array>
#include <utility>

#include "../BoundaryModels/BounceBack.hh"
#include "../Collide.hh"
#include "../Data.hh"
#include "../Forces/ChemicalForce.hh"
#include "../Geometry.hh"
#include "../GradientStencils/GradientStencils.hh"
#include "../Parameters.hh"
#include "FlowField.hh"
#include "ModelBase.hh"

// Binary.hh: Contains the details of the LBM model to solve an equation for phase separation. Each
// Model is given a "TTraits" class that contains stencil, data, force and boundary information

template <class TLattice>
using DefaultTraitBinary = typename DefaultTrait<TLattice, 2>::template SetBoundary<BounceBack>::
    template SetProcessor<GradientsMultiStencil<OrderParameter<>, CentralXYZ, LaplacianCentral>>::template AddProcessor<
        std::tuple<ChemicalPotentialCalculatorBinary, CubicWetting>>;

template <class TLattice, class TTraits = DefaultTraitBinary<TLattice>>
class Binary : public CollisionBase<TLattice, typename TTraits::Stencil>,
               public ModelBase<TLattice, TTraits> {  // Inherit from base class to avoid repetition of common
    // calculations

    using Stencil = typename TTraits::Stencil;
    static constexpr int mNDIM = TLattice::NDIM;

   public:
    inline void setTau1(double val) { mTau1 = val; }
    inline void setTau2(double val) { mTau2 = val; }
    inline void setA(double val) { mA = val; }

    inline void collide() override;  // Collision step

    inline void initialise() override;  // Initialisation step

    inline void computeMomenta() override;  // Momenta (density, velocity) calculation

    inline double computeEquilibrium(int k, int idx) override;  // Calculate equilibrium in direction idx

    /**
     * \brief Function to compute the fractional concentration of each fluid (between 0 and 1).
     * \param k The index of the node on the lattice.
     * \param iFluid The index of the fluid (0 or 1).
     */
    inline double computeConcentration(int k, int iFluid) override;

   private:
    inline double computeModelForce(int xyz, int k);  // Calculate forces specific to the model in direction xyz

    static constexpr double mTau = 1.0;                // TEMPORARY relaxation time
    static constexpr double mInverseTau = 1.0 / mTau;  // TEMPORARY inverse relaxation time

    enum { x = 0, y = 1, z = 2 };  // Indices corresponding to x, y, z directions

    double mGamma = 1;

    double mTau1 = 1;
    double mTau2 = 1;

    double mA;

    std::vector<double>& orderparameter = OrderParameter<>::get<TLattice>();  // Reference to vector of order parameters
    std::vector<double>& velocity = Velocity<>::get<TLattice, TLattice::NDIM>();  // Reference to vector of velocities
    std::vector<double>& itau = InverseTau<>::get<TLattice>();    // Reference to vector of inverse relaxation times
    std::vector<double>& pressure = Pressure<>::get<TLattice>();  // Reference to vector of pressures
};

template <class TLattice, class TTraits>
inline void Binary<TLattice, TTraits>::collide() {
#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

        if (this->isCollisionNode(k)) {
            double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

            double sum = 0;

            double equilibriums[Stencil::Q] = {};

            for (int idx = 1; idx < Stencil::Q; idx++) {
                equilibriums[idx] = computeEquilibrium(k, idx);
                sum += equilibriums[idx];
            }

            equilibriums[0] = orderparameter[k] - sum;

            this->collisionQ(equilibriums, old_distribution, mInverseTau, k);
        }
    }

    this->mData.communicateDistribution();
}

template <class TLattice, class TTraits>
inline void Binary<TLattice, TTraits>::initialise() {  // Initialise model
    this->initialiseProcessors();

    this->mData.generateNeighbors();  // Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this->mt_Forces, mTau1, mTau2);

#pragma omp parallel for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

        double* distribution = this->mDistribution.getDistributionPointer(k);
        double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

        ChemicalPotential<>::initialise<TLattice>(0, k);

        OrderParameter<>::initialise<TLattice>(1.0, k);

        // TODO: Decide what this should be
        // InverseTau = c_1 / tau_1 + c_2 / tau_2; c_1 = 0.5 * (rho + phi), c_2 = 0.5 * (rho - phi)
        double density = Density<>::get<TLattice>(k);
        double orderParameter = OrderParameter<>::get<TLattice>(k);
        InverseTau<>::initialise<TLattice>(
            0.5 * ((density + orderParameter) / mTau1 + (density - orderParameter) / mTau2), k);
        // InverseTau<>::initialise<TLattice>( 1.0 / (0.5 * (1.0 + orderParameter) * (mTau1 - mTau2) + mTau2), k);

        Pressure<>::initialise<TLattice>(density / 3., k);

        double equilibriumsum = 0;

        for (int idx = Stencil::Q - 1; idx >= 0; idx--) {
            double equilibrium;

            if (idx > 0)
                equilibrium = computeEquilibrium(k, idx);
            else
                equilibrium = orderparameter[k] - equilibriumsum;

            distribution[idx] = equilibrium;  // Set distributions to equillibrium
            old_distribution[idx] = equilibrium;

            equilibriumsum += equilibrium;
        }
    }

    this->mData.communicate(BoundaryLabels<TTraits::Lattice::NDIM>::template getInstance<TLattice>());
}

template <class TLattice, class TTraits>
inline void Binary<TLattice, TTraits>::computeMomenta() {  // Calculate order parameter

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // Loop over k

        if (this->isCollisionNode(k)) {
            double* distribution = this->mDistribution.getDistributionPointer(k);

            orderparameter[k] = this->computeDensity(distribution, k);

            double density = Density<>::get<TLattice>(k);

            // TODO: Decide what this should be
            itau[k] = 0.5 * ((density + orderparameter[k]) / mTau1 + (density - orderparameter[k]) / mTau2);
            // itau[k] = 1.0 / (0.5 * (1.0 + orderparameter[k]) * (mTau1)-0.5 * (-1.0 + orderparameter[k]) * mTau2);

            pressure[k] = density / 3 + mA * (-0.5 * pow(orderparameter[k], 2) + 0.75 * pow(orderparameter[k], 4));
        }
    }

    this->mData.communicate(OrderParameter<>::getInstance<TLattice>());
}

template <class TLattice, class TTraits>
inline double Binary<TLattice, TTraits>::computeEquilibrium(int k, int idx) {
    double velocityFactor = CollisionBase<TLattice, Stencil>::computeVelocityFactor(&velocity[k * Stencil::D], idx);
    return Stencil::Weights[idx] *
           (ChemicalPotential<>::get<TLattice>(k) * mGamma / Stencil::Cs2 + orderparameter[k] * velocityFactor);
}

template <class TLattice, class TTraits>
inline double Binary<TLattice, TTraits>::computeConcentration(int k, int iFluid) {
    if (iFluid == 0) {
        return 0.5 * (1.0 + orderparameter[k]);
    } else if (iFluid == 1) {
        return 0.5 * (1.0 - orderparameter[k]);
    } else {
        throw std::invalid_argument("Invalid fluid number. Must be 0 or 1.");
    }
}

// FlowField.hh: Contains the details of the LBM model to solve the Navier-Stokes and continuity equation. Each
// Model is given a "TTraits" class that contains stencil, data, force and boundary information

template <class TLattice>
using DefaultTraitFlowFieldBinary =
    typename DefaultTrait<TLattice,
                          2>::template SetBoundary<BounceBack>::template SetForce<ChemicalForceBinary<Guo, Gradient>>;

template <class TLattice, class TTraits = DefaultTraitFlowFieldBinary<TLattice>>
class FlowFieldBinary : public FlowField<TLattice, TTraits> {  // Inherit from base class to avoid repetition of common
    // calculations

    using Stencil = typename TTraits::Stencil;
    static constexpr int mNDIM = TLattice::NDIM;

   public:
    inline void setTau1(double val) { mTau1 = val; }
    inline void setTau2(double val) { mTau2 = val; }

    inline virtual void collide() override;  // Collision step

    inline virtual void initialise() override;  // Initialisation step

    inline double computeEquilibrium(int k, int idx) override;  // Calculate equilibrium in direction idx

    /**
     * \brief Function to compute the fractional concentration of each fluid (between 0 and 1).
     * \param k The index of the node on the lattice.
     * \param iFluid The index of the fluid (0 or 1).
     */
    inline double computeConcentration(int k, int iFluid) override;

   private:
    double mTau1 = 1;
    double mTau2 = 1;

    enum { x = 0, y = 1, z = 2 };
};

template <class TLattice, class TTraits>
inline double FlowFieldBinary<TLattice, TTraits>::computeEquilibrium(int k, int idx) {
    double density = this->density[k];
    double order_parameter = OrderParameter<>::get<TLattice>(k);
    double chemical_potential = ChemicalPotential<>::get<TLattice>(k);
    double gamma = CollisionBase<TLattice, Stencil>::computeGamma(&(this->velocity[k * Stencil::D]), idx);
    return density * gamma + Stencil::Weights[idx] * order_parameter * chemical_potential /
                                 Stencil::Cs2;  // Equilibrium is density times gamma in this case
}

template <class TLattice, class TTraits>
inline void FlowFieldBinary<TLattice, TTraits>::collide() {  // Collision step

#pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

        if (this->isCollisionNode(k)) {
            double* old_distribution = this->mDistribution.getDistributionOldPointer(k);
            double equilibriumsum = 0;

            double equilibriums[Stencil::Q];

            for (int idx = 1; idx < Stencil::Q; idx++) {
                equilibriums[idx] = computeEquilibrium(k, idx);
                equilibriumsum += equilibriums[idx];
            }

            equilibriums[0] = this->density[k] - equilibriumsum;

            this->collisionQ(equilibriums, old_distribution, InverseTau<>::get<TLattice>(k), k);
        }
    }

    FlowField<TLattice, TTraits>::mData.communicateDistribution();
}

template <class TLattice, class TTraits>
inline void FlowFieldBinary<TLattice, TTraits>::initialise() {  // Initialise model
    this->initialiseProcessors();

    this->mData.generateNeighbors();  // Fill array of neighbor values (See Data.hh)
    TTraits::template CollisionModel<Stencil>::template initialise<TLattice>(this->mt_Forces, mTau1, mTau2);

#pragma omp parallel for schedule(guided)
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {  // loop over k

        double* distribution = this->mDistribution.getDistributionPointer(k);
        double* old_distribution = this->mDistribution.getDistributionOldPointer(k);

        Density<>::initialise<TLattice>(1.0, k);  // Set density to 1 initially (This will change)
        Velocity<>::initialise<TLattice, TLattice::NDIM>(0.0, k, x);
        if constexpr (TLattice::NDIM >= 2) Velocity<>::initialise<TLattice, TLattice::NDIM>(0.0, k, y);
        if constexpr (TLattice::NDIM == 3) Velocity<>::initialise<TLattice, TLattice::NDIM>(0.0, k, z);

        double equilibriumsum = 0;

        for (int idx = Stencil::Q - 1; idx >= 0; idx--) {
            double equilibrium;

            if (idx > 0)
                equilibrium = computeEquilibrium(k, idx);

            else
                equilibrium = FlowField<TLattice, TTraits>::density[k] - equilibriumsum;

            distribution[idx] = equilibrium;  // Set distributions to equillibrium
            old_distribution[idx] = equilibrium;

            equilibriumsum += equilibrium;
        }
    }
}

template <class TLattice, class TTraits>
inline double FlowFieldBinary<TLattice, TTraits>::computeConcentration(int k, int iFluid) {
    double orderParameter = OrderParameter<>::get<TLattice>(k);
    if (iFluid == 0) {
        return 0.5 * (1.0 + orderParameter);
    } else if (iFluid == 1) {
        return 0.5 * (1.0 - orderParameter);
    } else {
        throw std::invalid_argument("Invalid fluid number. Must be 0 or 1.");
    }
}
