#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralQ : GradientBase<AllDirections> {

    template<class T_traits, class T_parameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class T_obj>
    using GradientType = Gradient<T_obj,T_obj::instances>;
    
};

template<class T_traits,class T_parameter>
inline double CentralQ::compute(const int direction, const int k, const int num) {

    using Lattice = typename T_traits::Lattice;
    using Stencil = typename T_traits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();
    
    return 0.5 * (T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) - T_parameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));
        
}