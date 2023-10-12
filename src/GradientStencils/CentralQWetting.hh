#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQWetting : WettingGradient<AllDirections> {

    template<class TTraits, class TParameter>
    inline double compute(const int direction, const int k, int num = 0);

    template<class TObj>
    using GradientType = Gradient<TObj,TObj::instances>;
    
};

template<class TTraits,class TParameter>
inline double CentralQWetting::compute(const int direction, const int k, const int num) {

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (Geometry<Lattice>::getBoundaryType(k) == 4) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction]) == 1)) {

        double csolid = TParameter::template get<Lattice>(k, num);//TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalq), num);

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]) == 1)) return 0;
        return 0.5 * ((csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))) - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]) == 1)) {

        double csolid = TParameter::template get<Lattice>(k, num);//TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalq), num);

        return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) - (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))));

    }
    else {

        return 0.5 * (TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
        
}