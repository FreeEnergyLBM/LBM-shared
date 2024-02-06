#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedQBounceBack : GradientBase<AllDirections> {

    template<class TTraits, class TParameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class TObj>
    using GradientType = GradientMixed<TObj,TObj::instances>;
    
};

template<class TTraits, class TParameter>
inline double MixedQBounceBack::compute(const int direction, const int k, int num){

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    /*
    return 0.25 * (-TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction], num)
                   + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) 
                   - 3 * TParameter::template get<Lattice>(k, num) 
                   - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));
    */
    
    const static auto& param = TParameter::template get<Lattice>();
    /*
    if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])==1)) {

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)){
            return 0;
        }

        return 0.25 * (2 *  TParameter::template get<Lattice>(k, num)
                       - 2 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction])==1)) {

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)){
            return 0.25 * (4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q+  direction], num) 
                       - 4 * TParameter::template get<Lattice>(k, num));
        }

        return 0.25 * (4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q+  direction], num) 
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    
    else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)) {

        return 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction], num)
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                       - 4 * TParameter::template get<Lattice>(k, num));

    }
    else {

        return 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q+direction], num)
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    */
    if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])!=1)
                && (Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])!=1)) {

        return 0.25 * (- param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q+direction]*TParameter::instances + num]
                       + 5 * param[data.getNeighbors()[k * Stencil::Q + direction]*TParameter::instances + num]
                       - 3 * param[k*TParameter::instances + num]
                       - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]);

    }
    else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])==1)) {

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)){
            return 0;
        }

        return 0.25 * (2 *  param[k*TParameter::instances + num]
                       - 2 * param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]);

    }
    else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction])==1)) {

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)){
            return 0.25 * (4 * param[data.getNeighbors()[k * Stencil::Q+  direction]*TParameter::instances + num] 
                       - 4 * param[k*TParameter::instances + num]);
        }

        return 0.25 * (4 * param[data.getNeighbors()[k * Stencil::Q+  direction]*TParameter::instances + num] 
                       - 3 * param[k*TParameter::instances + num]
                       - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]);

    }
    
    else {

        return 0.25 * (- param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction]*TParameter::instances + num]
                       + 5 * param[data.getNeighbors()[k * Stencil::Q + direction]*TParameter::instances + num]
                       - 4 * param[k*TParameter::instances + num]);

    }
    
}