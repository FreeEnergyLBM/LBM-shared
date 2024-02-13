#pragma once
#include "AddOnBase.hh"
#include <iostream>

template<class TParam>
class NoFluxSolid : public AddOnBase {
    public:

        NoFluxSolid() { this->setNodeID(5); }

        template<class TTraits>
        inline void compute(int k);

};

template<class TParam>
template<class TTraits>
inline void NoFluxSolid<TParam>::compute(int k){

    if (Geometry<typename TTraits::Lattice>::isBulkSolid(k)) return;

    //int normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection)->second;

    using Lattice = typename TTraits::Lattice;
    using DataType = Data_Base<Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    if (!this->apply<Lattice>(k)) return;

    /*
    bool cont = true;
    bool cont2 = true;

    for (int i : mNodeID){
        if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(data.getNeighbor(k, normalq), normalq)) == i) goto dontapply;
        else if(Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, normalq)) == i) goto applytwoneighbors;
    }

    cont = false;
    cont2 = false;

    dontapply:
        if (cont) return;

    applytwoneighbors:
    //#pragma omp critical
    //{
    //std::cout<<k<<" "<<normalq<<std::endl;
    //}
    if (cont2) {
        TParam::template get<Lattice>(k) = TParam::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, normalq), normalq));
        //count += 1.0;
    }
    else {
        TParam::template get<Lattice>(k) = TParam::template get<Lattice>(data.getNeighbor(k, normalq));
        //count += 1.0;
    }
    */

    double val = 0;
    double count = 0;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

        int kNeighbor1 = data.getNeighbor(k, idx);
        int kNeighbor2 = data.getNeighbor(kNeighbor1, idx);
        if (this->apply<Lattice>(kNeighbor2)) {
            continue;
        } else if (this->apply<Lattice>(kNeighbor1)) {
            val += TParam::template get<Lattice>(kNeighbor2);
        } else {
            val += TParam::template get<Lattice>(kNeighbor1);
        }
        count += 1.0;

    }

    if(count>0) TParam::template get<Lattice>(k) = val/count;

}
