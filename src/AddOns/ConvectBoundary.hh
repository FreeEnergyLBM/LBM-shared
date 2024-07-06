#pragma once
#include <math.h>

#include "../Parameters.hh"
#include "AddOnBase.hh"

template <class TParameter, class TParameterOld>
class ConvectParameterBoundary : public AddOnBase {
   public:
    ConvectParameterBoundary() { this->setNodeID(4); }

    template <class TTraits>
    inline void compute(int k);

    template <class TTraits>
    inline void communicate();

    double mVelocity = 0;
    int mCount = 1;
};

template <class TParameter, class TParameterOld>
template <class TTraits>
inline void ConvectParameterBoundary<TParameter, TParameterOld>::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    const int NDIM = Lattice::NDIM;
    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    if (!this->apply<Lattice>(k)) return;

    const std::array<int8_t, NDIM>& normal = BoundaryLabels<NDIM>::template get<Lattice>(k).NormalDirection;
    int normalq = Stencil::QMap.find(normal)->second;
    if (normalq > 0) {
        double normalvelocity = 0;
        double magnormal = 0;

        for (int xyz = 0; xyz < NDIM; xyz++) {
            normalvelocity += -normal[xyz] * Velocity<>::get<Lattice, NDIM>(
                                                 data.getNeighbor(data.getNeighbor(k, normalq), normalq), xyz);
            magnormal += pow(normal[xyz], 2);
        }

        magnormal = sqrt(magnormal);

        normalvelocity *= 1. / magnormal;
        normalvelocity = mVelocity / ((double)mCount);

        std::vector<double>& param = TParameter::template get<Lattice>();
        const std::vector<double>& paramOld = TParameterOld::template get<Lattice>();
        param[k] = (paramOld[k] + normalvelocity * param[data.getNeighbor(k, normalq)]) / (1 + normalvelocity);
        param[data.getNeighbor(k, Stencil::Opposites[normalq])] =
            (paramOld[data.getNeighbor(k, Stencil::Opposites[normalq])] + normalvelocity * param[k]) /
            (1 + normalvelocity);
        //if (TIME==40000) std::cout<<param[data.getNeighbor(k, normalq)]<<" "<<param[k]<<" "<<param[data.getNeighbor(k, Stencil::Opposites[normalq])]<<std::endl;

    }
}

template <class TParameter, class TParameterOld>
template <class TTraits>
inline void ConvectParameterBoundary<TParameter, TParameterOld>::communicate() {
    using Lattice = typename TTraits::Lattice;
    #pragma omp master
    {
        mVelocity = 0;
        mCount = 0;
    }
    #pragma omp barrier
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    #pragma omp parallel for schedule(guided)
    for (int k = 0; k < Lattice::N; k++) {

        if (!this->apply<Lattice>(k)) continue;

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
}
