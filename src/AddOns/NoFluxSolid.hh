#pragma once
#include<iostream>

template<class TParam>
class NoFluxSolid : public AddOnBase {
    
    public:

        template<class TTraits>
        inline void compute(int k);

        inline void setInterfaceID(int id) {mInterfaceID[0]=id;};

        inline void setInterfaceID(const std::vector<int>& id) {mInterfaceID=id;};

    private:

        std::vector<int> mInterfaceID = {1};

};

template<class TParam>
template<class TTraits>
inline void NoFluxSolid<TParam>::compute(int k){

    int normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection)->second;

    using DataType = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    for (int i : mInterfaceID){
        if(Geometry<typename TTraits::Lattice>::getBoundaryType(k) == i && Geometry<typename TTraits::Lattice>::getBoundaryType(data.getNeighbor(data.getNeighbor(k, normalq), normalq)) != i) goto runloop;
    }

    return;

    runloop:

        TParam::template get<typename TTraits::Lattice>(k) = TParam::template get<typename TTraits::Lattice>(data.getNeighbor(k, normalq));

}