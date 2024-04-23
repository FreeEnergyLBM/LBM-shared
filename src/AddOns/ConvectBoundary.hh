#pragma once
#include "../Parameters.hh"
#include "AddOnBase.hh"
#include <iostream>
#include <math.h>

template<class TParameter, class TParameterOld>
class ConvectParameterBoundary : public AddOnBase {
    public:

        ConvectParameterBoundary() { this->setNodeID(4); }

        template<class TTraits>
        inline void compute(int k);

};

template<class TParameter, class TParameterOld>
template<class TTraits>
inline void ConvectParameterBoundary<TParameter,TParameterOld>::compute(int k) {

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    if (!this->apply<Lattice>(k)) return;

    //const std::vector<int>& neighbors = DataType::getInstance().getNeighbors();
    const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;
    int normalq = Stencil::QMap.find(normal)->second;
    if(normalq>0){

        double normalvelocity = 0;
        double magnormal = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            normalvelocity += -normal[xyz]*Velocity<>::get<Lattice, Lattice::NDIM>(DataType::getInstance().getNeighbor(DataType::getInstance().getNeighbor(k,normalq),normalq),xyz);
            magnormal += pow(normal[xyz],2);
        }

        magnormal = sqrt(magnormal);

        normalvelocity *= 1./magnormal;
        for (int component = 0 ; component < TParameter::instances; component++){
            TParameter::template get<Lattice>(k,component) = (TParameterOld::template get<Lattice>(k,component)+normalvelocity*TParameter::template get<Lattice>(DataType::getInstance().getNeighbor(k,normalq),component))/(1+normalvelocity);

            TParameter::template get<Lattice>(data.getNeighbor(k, Stencil::Opposites[normalq]),component) = (TParameterOld::template get<Lattice>(data.getNeighbor(k, Stencil::Opposites[normalq]),component)+normalvelocity*TParameter::template get<Lattice>(k,component))/(1+normalvelocity);
        }

    }

}
