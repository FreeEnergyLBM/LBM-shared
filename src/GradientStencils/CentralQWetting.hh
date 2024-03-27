#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQWetting : GradientBase<Gradient,AllDirections> {

    template<class TTraits, class TParameter>
    inline double compute(const int direction, const int k, int num = 0);
    
};
/*
template<class TTraits,class TParameter>
inline double CentralQWetting::compute(const int direction, const int k, const int num) {

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (this->isBoundary<Lattice>(k) == 4) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]) == 1)) {

        double csolid = TParameter::template get<Lattice>(k, num);//TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalq), num);

        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]) == 1)) return 0;
        return 0.5 * ((csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))) - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]) == 1)) {

        double csolid = TParameter::template get<Lattice>(k, num);//TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalq), num);

        return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) - (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))));

    }
    else {

        return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
        
}
*/

template<class TTraits, class TParameter>
inline double CentralQWetting::compute(const int direction, const int k, int num){
        
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    //if (this->isBoundary<Lattice>(k) == 4) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction]))) {

        const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, direction)).NormalDirection)->second;

        double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalq), num);

        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {
            const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;
            double csolidbackward = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalqbackward), num);
            return 0.5 * ((csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))) - (csolidbackward - 0.5 * this->mPrefactor * (csolidbackward - pow(csolidbackward, 2))));
        }
        return 0.5 * ((csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))) - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {

        const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;

        double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalq), num);

        return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) - (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))));

    }
    else {

        return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
        
}