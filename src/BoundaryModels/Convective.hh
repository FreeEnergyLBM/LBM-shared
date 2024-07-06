#pragma once
#include <iostream>

#include "../Parameters.hh"
#include "BoundaryBase.hh"

class Convective : public BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);

    template <class TTraits>
    inline void communicate(){};

    template <class TTraits, class TDistributionType>
    inline void communicate(TDistributionType& mDistribution);

    template <class TTraits>
    inline void communicateProcessor();

    template <class TTraits>
    inline void runProcessor(int k);
    double mVelocity = 0;
    int mCount = 0;
};

template <class TTraits>
inline void Convective::communicateProcessor() {
#pragma omp single
    {
        mVelocity = 0;
        mCount = 0;
    }
}

template <class TTraits>
inline void Convective::runProcessor(int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    using DataType = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    const int& normalq =
        Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;
    const std::array<int8_t, TTraits::Lattice::NDIM>& normal =
        BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;

    double normalvelocity = 0;
    double magnormal = 0;

    for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
        normalvelocity += -normal[xyz] * (Velocity<>::get<Lattice, Lattice::NDIM>(
                                             data.getNeighbor(data.getNeighbor(k, normalq), normalq), xyz));
        magnormal += pow(normal[xyz], 2);
    }

    magnormal = sqrt(magnormal);

    normalvelocity *= 1. / magnormal;
#pragma omp critical
    {
        mVelocity += normalvelocity;
        mCount += 1;
    }
}

template <class TTraits, class TDistributionType>
inline void Convective::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    const int& normalq =
        Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;
    const std::array<int8_t, TTraits::Lattice::NDIM>& normal =
        BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        double normdotci = 0;
        for (int xyz = 0; xyz < Lattice::NDIM; xyz++) {
            normdotci += normal[xyz] * Stencil::Ci_xyz(xyz)[idx];
        }
        if (normdotci <= 0) continue;

        double normalvelocity = 0;
        double magnormal = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            normalvelocity +=
                -normal[xyz] * Velocity<>::get<Lattice, Lattice::NDIM>(
                                   distribution.streamIndex(distribution.streamIndex(k, normalq), normalq), xyz);
            magnormal += pow(normal[xyz], 2);
        }

        magnormal = sqrt(magnormal);

        normalvelocity *= 1. / magnormal;
        normalvelocity = mVelocity / ((double)mCount);

        distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] =
            (distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[idx] +
             normalvelocity * distribution.getDistributionPointer(
                                  distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[idx]) /
            (1 + normalvelocity);
    }
}

template <class TTraits, class TDistributionType>
inline void Convective::communicate(TDistributionType& distribution) {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
}

template <class TForceTuple>
class Convective2 : public BoundaryBase {
   public:
    using ForceTuple = TForceTuple;
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);

    template <class TTraits>
    inline void communicate(){};

    inline void setForceTuple(const ForceTuple& tup) { mt_Forces = tup; }

    template <class TTraits, class TDistributionType>
    inline void communicate(TDistributionType& mDistribution);

    template <class TTraits>
    inline void communicateProcessor();

    template <class TTraits>
    inline void runProcessor(int k);

    inline void setVelocityCalculator(double (*v)(const double* distribution, ForceTuple& forcetuple,
                                                  const double& density, int xyz, int k)) {
        mVelocityCalculator = v;
    }

   private:
    TForceTuple mt_Forces;

    static double defaultVelocityCalculator(const double* distribution, ForceTuple& forcetuple, const double& density,
                                            int xyz, int k) {
        return 0;
    }

    double (*mVelocityCalculator)(const double* distribution, ForceTuple& forcetuple, const double& density, int xyz,
                                  int k) = &defaultVelocityCalculator;

    double mVelocity = 0;
    int mCount = 0;
};

template <class TForceTuple>
template <class TTraits>
inline void Convective2<TForceTuple>::runProcessor(int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    using DataType = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    const int& normalq =
        Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;
    const std::array<int8_t, TTraits::Lattice::NDIM>& normal =
        BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;

    double normalvelocity = 0;
    double magnormal = 0;

    for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
        normalvelocity += -normal[xyz] * (Velocity<>::get<Lattice, Lattice::NDIM>(
                                             data.getNeighbor(data.getNeighbor(k, normalq), normalq), xyz));
        magnormal += pow(normal[xyz], 2);
    }

    magnormal = sqrt(magnormal);

    normalvelocity *= 1. / magnormal;
#pragma omp critical
    {
        mVelocity += normalvelocity;
        mCount += 1;
    }
}

template <class TForceTuple>
template <class TTraits, class TDistributionType>
inline void Convective2<TForceTuple>::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    const int& normalq =
        Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;
    const std::array<int8_t, TTraits::Lattice::NDIM>& normal =
        BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;

    for (int idx = 0; idx < Stencil::Q; idx++) {
        double normdotci = 0;
        for (int xyz = 0; xyz < Lattice::NDIM; xyz++) {
            normdotci += normal[xyz] * Stencil::Ci_xyz(xyz)[idx];
        }
        if (normdotci <= 0) continue;
        // std::cout<<idx<<std::endl;

        double civelocity = 0;
        double normalvelocity = 0;
        double magnormal = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            civelocity +=
                Stencil::Ci_xyz(xyz)[idx] *
                Density<>::get<Lattice>(distribution.streamIndex(
                    distribution.streamIndex(distribution.streamIndex(k, normalq), normalq), normalq)) *
                (3 * mVelocityCalculator(distribution.getDistributionPointer(distribution.streamIndex(k, normalq)),
                                         mt_Forces, Density<>::get<Lattice>(distribution.streamIndex(k, normalq)), xyz,
                                         distribution.streamIndex(k, normalq)) -
                 4 * mVelocityCalculator(distribution.getDistributionPointer(
                                             distribution.streamIndex(distribution.streamIndex(k, normalq), normalq)),
                                         mt_Forces,
                                         Density<>::get<Lattice>(
                                             distribution.streamIndex(distribution.streamIndex(k, normalq), normalq)),
                                         xyz, distribution.streamIndex(distribution.streamIndex(k, normalq), normalq)) +
                 mVelocityCalculator(
                     distribution.getDistributionPointer(distribution.streamIndex(
                         distribution.streamIndex(distribution.streamIndex(k, normalq), normalq), normalq)),
                     mt_Forces,
                     Density<>::get<Lattice>(distribution.streamIndex(
                         distribution.streamIndex(distribution.streamIndex(k, normalq), normalq), normalq)),
                     xyz,
                     distribution.streamIndex(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq),
                                              normalq)));
            magnormal += pow(normal[xyz], 2);
            normalvelocity +=
                -normal[xyz] *
                mVelocityCalculator(
                    distribution.getDistributionPointer(
                        distribution.streamIndex(distribution.streamIndex(k, normalq), normalq)),
                    mt_Forces,
                    Density<>::get<Lattice>(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq)),
                    xyz, distribution.streamIndex(distribution.streamIndex(k, normalq), normalq));
        }

        magnormal = sqrt(magnormal);

        normalvelocity *= 1. / magnormal;

        distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] =
            (distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[idx] +
             normalvelocity * distribution.getDistributionPointer(
                                  distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[idx]) /
            (1.0 + normalvelocity);
    }
}

template <class TForceTuple>
template <class TTraits, class TDistributionType>
inline void Convective2<TForceTuple>::communicate(TDistributionType& distribution) {
    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
}

template <class TForceTuple>
template <class TTraits>
inline void Convective2<TForceTuple>::communicateProcessor() {
#pragma omp single
    {
        mVelocity = 0;
        mCount = 0;
    }
}
