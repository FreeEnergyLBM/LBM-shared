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

    if (Geometry<typename TTraits::Lattice>::isBulkSolid(k)) return;

    //int normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection)->second;

    using DataType = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    for (int i : mInterfaceID){
        if(Geometry<typename TTraits::Lattice>::getBoundaryType(k) == i) goto runloop;// && Geometry<typename TTraits::Lattice>::getBoundaryType(data.getNeighbor(data.getNeighbor(k, normalq), normalq)) != i) goto runloop;
    }

    return;

    runloop:
        /*
        bool cont = true;
        bool cont2 = true;

        for (int i : mInterfaceID){
            if (Geometry<typename TTraits::Lattice>::getBoundaryType(data.getNeighbor(data.getNeighbor(k, normalq), normalq)) == i) goto dontapply;
            else if(Geometry<typename TTraits::Lattice>::getBoundaryType(data.getNeighbor(k, normalq)) == i) goto applytwoneighbors;
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
            TParam::template get<typename TTraits::Lattice>(k) = TParam::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, normalq), normalq));
            //count += 1.0;
        }
        else {
            TParam::template get<typename TTraits::Lattice>(k) = TParam::template get<typename TTraits::Lattice>(data.getNeighbor(k, normalq));
            //count += 1.0;
        }
        */
        
        double val = 0;
        double count = 0;

        for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

            bool cont = true;
            bool cont2 = true;

            for (int i : mInterfaceID){
                if (Geometry<typename TTraits::Lattice>::getBoundaryType(data.getNeighbor(data.getNeighbor(k, idx), idx)) == i) goto dontapply;
                else if(Geometry<typename TTraits::Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == i) goto applytwoneighbors;
            }

            cont = false;
            cont2 = false;

            dontapply:
                if (cont) continue;

            applytwoneighbors:

            if (cont2) {
                val += TParam::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx));
                count += 1.0;
            }
            else {
                val += TParam::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx));
                count += 1.0;
            }
            
            
        }

        if(count>0) TParam::template get<typename TTraits::Lattice>(k) = val/count;
        

}