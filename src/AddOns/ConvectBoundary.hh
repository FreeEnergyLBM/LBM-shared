#pragma once
#include "../Parameters.hh"
#include "AddOnBase.hh"
#include<iostream>
#include<math.h>

template<class TParameter, class TParameterOld>
class ConvectParameterBoundary : public AddOnBase {
    public:

        template<class TTraits>
        inline void compute(int k);

        inline void setInterfaceID(int id) {mInterfaceID[0]=id;};

        inline void setInterfaceID(const std::vector<int>& id) {mInterfaceID=id;};

    private:

        double mInterfaceVal;
        std::vector<int> mInterfaceID = {4};

};

template<class TParameter, class TParameterOld>
template<class TTraits>
inline void ConvectParameterBoundary<TParameter,TParameterOld>::compute(int k) {

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    for (int i : mInterfaceID){
        if(Geometry<typename TTraits::Lattice>::getBoundaryType(k) == i) goto runloop;
    }

    return;

    runloop:

        const std::vector<int>& neighbors = DataType::getInstance().getNeighbors();
        const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;
        int normalq = Stencil::QMap.find(normal)->second;
        if(normalq>0){

            double normalvelocity = 0;
            double magnormal = 0;
            for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
                normalvelocity += normal[xyz]*Velocity<>::get<Lattice, Lattice::m_NDIM>(k,xyz);
                magnormal += pow(normal[xyz],2);
            }

            magnormal = sqrt(magnormal);

            normalvelocity *= 1./magnormal;

            TParameter::template get<Lattice>(k) = (TParameterOld::template get<Lattice>(k)+0.5*normalvelocity*TParameter::template get<Lattice>(DataType::getInstance().getNeighbors(k,normalq)))/(1+0.5*normalvelocity);

            TParameter::template get<Lattice>(data.getNeighbor(k, Stencil::Opposites[normalq])) = (TParameterOld::template get<Lattice>(data.getNeighbor(k, Stencil::Opposites[normalq]))+0.5*normalvelocity*TParameter::template get<Lattice>(k))/(1+0.5*normalvelocity);

        }
    
}