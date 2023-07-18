#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQ : GradientBase<AllDirections> {

    template<class T_traits, class T_parameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class T_obj>
    using GradientType = Gradient<T_obj,T_obj::instances>;
    
};

template<class T_traits, class T_parameter>
inline double CentralQNoSolid::compute(const int direction, const int k, int num){
        
    using Lattice = typename T_traits::Lattice;
    using Stencil = typename T_traits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    if ((Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + direction]))) {

        return 0.5 * (T_parameter::template get<Lattice>(k, num) - T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    else if ((Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {

        return 0.5 * (T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q+direction], num)-T_parameter::template get<Lattice>(k, num));

    }
    else if ((Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + direction])) && (Geometry<Lattice>::isSolid(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]))) {

        return 0;

    }
    else {

        return 0.5 * (T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)-T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
        
}