#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQMirrorSolid : GradientBase<AllDirections> {

    template<class TTraits, class TParameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class TObj>
    using GradientType = Gradient<TObj,TObj::instances>;
    
};

template<class TTraits, class TParameter>
inline double CentralQMirrorSolid::compute(const int direction, const int k, int num){
        
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (Geometry<Lattice>::getBoundaryType(k) == 4) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();
    const static auto& param = TParameter::template get<Lattice>();
    /*
    if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction]) == 1)) {

        const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, direction)).NormalDirection)->second;

        double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalq), num);

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]) == 1)) {
            const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;
            double csolidbackward = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalqbackward), num);
            return 0.5 * (csolid - csolidbackward);
        }
        return 0.5 * (csolid - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]) == 1)) {

        const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;

        double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalq), num);

        return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) - csolid);

    }
    else {

        return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    */
    if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction]) == 1)) {

        const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, direction)).NormalDirection)->second;

        double csolid = param[data.getNeighbor(data.getNeighbor(k, direction), normalq)*TParameter::instances + num];

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]) == 1)) {
            const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;
            double csolidbackward = param[data.getNeighbor(data.getNeighbor(k, direction), normalqbackward)*TParameter::instances + num];
            return 0.5 * (csolid - csolidbackward);
        }
        return 0.5 * (csolid - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]);

    }
    else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]) == 1)) {

        const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;

        double csolid = param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalq)*TParameter::instances + num];

        return 0.5 * (param[data.getNeighbors()[k * Stencil::Q + direction]*TParameter::instances + num] - csolid);

    }
    else {

        return 0.5 * (param[data.getNeighbors()[k * Stencil::Q + direction]*TParameter::instances + num] - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]);

    }
}