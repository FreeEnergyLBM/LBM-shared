#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedQNoSolid : GradientBase<AllDirections> {

    template<class TTraits, class TParameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class TObj>
    using GradientType = Gradient<TObj,TObj::instances>;
    
};

template<class TTraits, class TParameter>
inline double MixedQNoSolid::compute(const int direction, const int k, int num){
        
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    if ((Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + direction]))) {

        return 0.25 * (2 *  TParameter::template get<Lattice>(k, num)
                       - 2 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    else if ((Geometry<Lattice>::isSolid(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]]))) {

        return 0.25 * (4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q+  direction], num) 
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    else if ((Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {

        return 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction], num)
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                       - 4 * TParameter::template get<Lattice>(k, num));

    }
    else if ((Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + direction]))
              && (Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {

        return 0;

    }
    else {

        return 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q+direction], num)
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
        
}