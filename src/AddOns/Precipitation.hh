#pragma once
#include <iostream>
#include <random>

#include "AddOnBase.hh"

template <class TParameter>
class Nucleation : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

    inline void setNucleationParameters(double criticalProbability, double criticalConcentration) {
        mNucleationActive = true;
        mCriticalNucleationProbability = criticalProbability;
        mCriticalSoluteConcentration = criticalConcentration;
    }

    inline void setNucleationSize(double size) { mNucleationSize = size; }

   private:
    template <class TTraits>
    double calculateNucleationProbability(int k, double randomValue);

    template <class TTraits>
    void performNucleation(int k, double nucleationSize);

    bool mNucleationActive = false;
    double mCriticalNucleationProbability = 0.8;
    double mCriticalSoluteConcentration = 0.6;
    double mNucleationSize = 1.0;  // Meaning one lattice unit in each direction (+/-)
};

template <class TParameter>
template <class TTraits>
inline void Nucleation<TParameter>::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    // using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    if (mNucleationActive) {
        // Generate random number with uniform distribution
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 0.5);
        double randomValue = dis(gen);

        double nucleationProbability = calculateNucleationProbability<TTraits>(k, randomValue);

        NucleationProbability<>::template get<Lattice>(k) = nucleationProbability;

        if (nucleationProbability > mCriticalNucleationProbability) {
            std::cout << "Nucleation at node " << k << " with probability " << nucleationProbability << std::endl;
            performNucleation<TTraits>(k, mNucleationSize);
        }
    }
}

template <class TParameter>
template <class TTraits>
double Nucleation<TParameter>::calculateNucleationProbability(int k, double randomValue) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    using data = Data_Base<Lattice, Stencil>;

    double nucleationProbability = randomValue;
    double orderParam = OrderParameter<>::template get<Lattice>(k);

    // We should have a higher probability of nucleation if being super-saturated
    // In the If-condition, we need to normlise the mCriticalSoluteConcentration with the order parameter
    if (TParameter::template get<Lattice>(k) >= mCriticalSoluteConcentration * orderParam) {
        double overSaturation = TParameter::template get<Lattice>(k) - mCriticalSoluteConcentration * orderParam;
        nucleationProbability += (0.1 + overSaturation);
    }

    if (orderParam <= 0.0) nucleationProbability *= 0.0;  // No nucleation in the gas phase

    // We should have a higher probability of nucleation if being adjacent to a solid boundary/salt crystal
    for (int idx = 1; idx < Stencil::Q; idx++) {
        int neighbor = data::getInstance().getNeighbors()[k * Stencil::Q + idx];
        if (Geometry<Lattice>::getBoundaryType(neighbor) == 1) nucleationProbability += 0.05;
    }

    // We should be away for already nucleated nodes
    auto [x, y, z] = computeXYZ<Lattice>(k);
    int investigationSize = 2;
    for (int i = -investigationSize; i <= investigationSize; i++) {
        for (int j = -investigationSize; j <= investigationSize; j++) {
            for (int l = -investigationSize; l <= investigationSize; l++) {
                int k_neighbor = computeKFromGlobal<Lattice>(x + i, y + j, z + l);
                if (NucleationProbability<>::template get<Lattice>(k_neighbor) == -1.0) {
                    nucleationProbability *= 0.0;
                    break;
                }
            }
        }
    }

    return nucleationProbability;
}

template <class TParameter>
template <class TTraits>
void Nucleation<TParameter>::performNucleation(int k, double nucleationSize) {
    using Lattice = typename TTraits::Lattice;

    // Modify the label of the fluid node to solid node in a square (2D)/cuboid (3D)of size 2*nucleationSize+1
    auto [x, y, z] = computeXYZ<Lattice>(k);
    for (int i = -nucleationSize; i <= nucleationSize; i++) {
        for (int j = -nucleationSize; j <= nucleationSize; j++) {
            for (int l = -nucleationSize; l <= nucleationSize; l++) {
                int k_neighbor = computeKFromGlobal<Lattice>(x + i, y + j, z + l);
                if (Geometry<Lattice>::getBoundaryType(k_neighbor) == 0) {
                    Geometry<Lattice>::modifyBoundaryLabel(k_neighbor, 1);
                    NucleationProbability<>::template get<Lattice>(k_neighbor) = -1.0;
                }
            }
        }
    }
}

class Growth : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);
    inline void setCriticalSoluteConcentration(double value) {
        mGrowthActive = true;
        mCriticalSoluteConcentration = value;
    }

   private:
    bool mGrowthActive = false;
    double mCriticalSoluteConcentration = 0.6;
};

template <class TTraits>
inline void Growth::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    using data = Data_Base<Lattice, Stencil>;

    if (!this->apply<Lattice>(k)) return;

    double orderParam = OrderParameter<>::template get<Lattice>(k);
    // mGrowthActive should be enabled and the order parameter should be >=0.8 (i.e. liquid phase) for the function
    // to consider the growth of the solid phase
    if (!mGrowthActive || orderParam <= 0.8) return;

    bool shouldModifyLabel = false;
    for (int idx = 1; idx < Stencil::Q; idx++) {
        int neighbor = data::getInstance().getNeighbors()[k * Stencil::Q + idx];
        if (NucleationProbability<>::template get<Lattice>(neighbor) == -1.0) {
            if (Solute<>::template get<Lattice>(neighbor) >= mCriticalSoluteConcentration) {
                std::cout << "Growth at node " << k << std::endl;
                shouldModifyLabel = true;
                break;
            }
        }
    }

    if (shouldModifyLabel) {
        Geometry<Lattice>::modifyBoundaryLabel(k, 1);
        NucleationProbability<>::template get<Lattice>(k) = -1.0;
    }
}

// The main purpose of this function is to modify label of the individual fluid nodes that become isolated at all
// directions due to precipitation to solid node
class InactiveIsolatedFluidNodes : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);
};

template <class TTraits>
inline void InactiveIsolatedFluidNodes::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    using data = Data_Base<Lattice, Stencil>;

    if (!this->apply<Lattice>(k)) return;

    int count = 0;
    for (int idx = 1; idx < Stencil::Q; idx++) {
        int neighbor = data::getInstance().getNeighbors()[k * Stencil::Q + idx];
        if (Geometry<Lattice>::getBoundaryType(neighbor) == 1) count++;
    }

    if (count >= Stencil::Q - 1) Geometry<Lattice>::modifyBoundaryLabel(k, 1);
}

/// ********** Old Precipitation model ********** ///

template <class TParameter>
class Precipitation : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

    inline void setPrecipitationThreshold(double value) {
        mPrecipitationActive = true;
        mPrecipitationThreshold = value;
    }

   private:
    template <class TTraits>
    double modifyPrecipitationProbability(int k, double randomValue);

    double mCriticalPreceipitationProbability = 0.8;
    double mPrecipitationThreshold = 1.0;
    bool mPrecipitationActive = false;
};

template <class TParameter>
template <class TTraits>
inline void Precipitation<TParameter>::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    // using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    // This is the simplest precipitation model, where precipitation occurs if the order parameter is above a certain
    // if (mPrecipitationActive && TParameter::template get<Lattice>(k) >= mPrecipitationThreshold) {
    //     Geometry<Lattice>::modifyBoundaryLabel(k, 1);
    // }

    if (mPrecipitationActive) {
        // Generate random number with normal distribution
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> dis(0.2, 0.1);  // mean 0.2, standard deviation 0.1
        double randomValue = dis(gen);

        double precipitationProbability = modifyPrecipitationProbability<TTraits>(k, randomValue);

        NucleationProbability<>::template get<Lattice>(k) = precipitationProbability;

        if (precipitationProbability > mCriticalPreceipitationProbability) {
            std::cout << "Precipitation at node " << k << " with probability " << precipitationProbability << std::endl;
            Geometry<Lattice>::modifyBoundaryLabel(k, 1);
        }
    }
}

template <class TParameter>
template <class TTraits>
double Precipitation<TParameter>::modifyPrecipitationProbability(int k, double randomValue) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    using data = Data_Base<Lattice, Stencil>;

    double precipitationProbability = randomValue;

    // We should have a higher probability of precipitation if being super-saturated
    if (TParameter::template get<Lattice>(k) >= mPrecipitationThreshold) {
        double overSaturation = TParameter::template get<Lattice>(k) - mPrecipitationThreshold;
        precipitationProbability += (0.2 + overSaturation);
    }

    // We should have a higher probability of precipitation if being closer to the liquid-vapour interface, i.e.
    // positive order parameters closer to zero
    if (OrderParameter<>::template get<Lattice>(k) < 0.0)
        precipitationProbability *= 0.0;  // No precipitation in the gas phase
    else {
        double orderParam = OrderParameter<>::template get<Lattice>(k);
        precipitationProbability += (orderParam < 0.8) ? 0.1 : -0.1;
    }

    // We should have a higher probability of precipitation if being adjacent to a solid boundary/salt crystal
    for (int idx = 1; idx < Stencil::Q; idx++) {
        int neighbor = data::getInstance().getNeighbors()[k * Stencil::Q + idx];
        if (Geometry<Lattice>::getBoundaryType(neighbor) == 1) precipitationProbability += 0.05;
    }

    return precipitationProbability;
}