#pragma once
#include "../Parameters.hh"
#include "AddOnBase.hh"
#include<iostream>
#include<math.h>

template<class TParameter>
class ConstantGradientBoundary : public AddOnBase {
    public:

        template<class TTraits>
        inline void compute(int k);

        inline void setInterfaceID(int id) {mInterfaceID[0]=id;};

        inline void setInterfaceID(const std::vector<int>& id) {mInterfaceID=id;};

    private:

        double mInterfaceVal;
        std::vector<int> mInterfaceID = {4};

};

template<class TParameter>
template<class TTraits>
inline void ConstantGradientBoundary<TParameter>::compute(int k) {

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    for (int i : mInterfaceID){
        if(Geometry<typename TTraits::Lattice>::getBoundaryType(k) == i) goto runloop;
    }

    return;

    runloop:

        //using data = Data_Base<Lattice, Stencil>;

        const std::vector<int>& neighbors = DataType::getInstance().getNeighbors();
        const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;
        int idx = Stencil::QMap.find(normal)->second;
        if(idx>0){

        double magnormal = 0;

        for (int xyz = 0; xyz < Lattice::NDIM; xyz++) {
            magnormal += normal[xyz]*normal[xyz];
        }

        magnormal = sqrt(magnormal);
        if (Geometry<Lattice>::isCorner(k)) {
            //TParameter::template get<Lattice>(neighbors[k * Stencil::Q+Stencil::Opposites[idx]]) = (4*TParameter::template get<Lattice>(k) - TParameter::template get<Lattice>(neighbors[k * Stencil::Q+idx])) / 3.0;
            
            for (int idx1 = 1; idx1 < Stencil::Q; idx1++) {
                double cidotnormal = 0;
                for (int xyz = 0; xyz < Lattice::NDIM; xyz++) {
                    cidotnormal += Stencil::Ci_xyz(xyz)[idx1]*normal[xyz];
                }
                if (cidotnormal > 0) {
                    TParameter::template get<Lattice>(k) = (4*TParameter::template get<Lattice>(neighbors[k * Stencil::Q+idx1]) - TParameter::template get<Lattice>(neighbors[neighbors[k * Stencil::Q+idx1] * Stencil::Q+idx1])) / 3.0;
                    TParameter::template get<Lattice>(neighbors[k * Stencil::Q+Stencil::Opposites[idx1]]) = (4*TParameter::template get<Lattice>(neighbors[k * Stencil::Q+idx1]) - TParameter::template get<Lattice>(neighbors[neighbors[k * Stencil::Q+idx1] * Stencil::Q+idx1])) / 3.0;
                }
            }
            /*
            if constexpr (Lattice::NDIM==)
            std::vector<int8_t> normal1 = {normal[0],(int8_t)0};
            int idx1 = Stencil::QMap.find(normal1)->second;
            TParameter::template get<Lattice>(neighbors[k * Stencil::Q+Stencil::Opposites[idx1]]) = (4*TParameter::template get<Lattice>(k) - TParameter::template get<Lattice>(neighbors[k * Stencil::Q+idx1])) / 3.0;

            std::vector<int8_t> normal2 = {(int8_t)0,normal[1]};
            int idx2 = Stencil::QMap.find(normal2)->second;
            TParameter::template get<Lattice>(neighbors[k * Stencil::Q+Stencil::Opposites[idx2]]) = (4*TParameter::template get<Lattice>(k) - TParameter::template get<Lattice>(neighbors[k * Stencil::Q+idx2])) / 3.0;
            */
        }
        else {
            TParameter::template get<Lattice>(k) = (4*TParameter::template get<Lattice>(neighbors[k * Stencil::Q+idx]) - TParameter::template get<Lattice>(neighbors[neighbors[k * Stencil::Q+idx] * Stencil::Q+idx])) / 3.0;

            TParameter::template get<Lattice>(data.getNeighbor(k, Stencil::Opposites[idx])) = (4*TParameter::template get<Lattice>(neighbors[k * Stencil::Q+idx]) - TParameter::template get<Lattice>(neighbors[neighbors[k * Stencil::Q+idx] * Stencil::Q+idx])) / 3.0;

        }

        }
    
}