#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQBounceBack : GradientBase<AllDirections> {

    template<class TTraits, class TParameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class TObj>
    using GradientType = Gradient<TObj,TObj::instances>;
    
};

template<class TTraits, class TParameter>
inline double CentralQBounceBack::compute(const int direction, const int k, int num){
        
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();
    const static auto& param = TParameter::template get<Lattice>();
    /*
    if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])==1)) {

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)) return 0;
        return 0.5 * (TParameter::template get<Lattice>(k, num) - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)) {

        return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q+direction], num)-TParameter::template get<Lattice>(k, num));

    }
    else {

        return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)-TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    */ 
    if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])==1)) {

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)) return 0;
        return 0.5 * (param[k*TParameter::instances + num] - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]);

    }
    else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)) {

        return 0.5 * (param[data.getNeighbors()[k * Stencil::Q+direction]*TParameter::instances + num]-param[k*TParameter::instances + num]);

    }
    else {

        return 0.5 * (param[data.getNeighbors()[k * Stencil::Q + direction]*TParameter::instances + num]-param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]);

    }
}