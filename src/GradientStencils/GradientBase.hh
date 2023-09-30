#pragma once
#include "../Service.hh"

template<class TDirections=Cartesian>
struct GradientBase{

    template<class TStencil>
    inline static constexpr int getNumberOfDirections(){
        if constexpr (std::is_same_v<TDirections,Cartesian>) return TStencil::D; 
        else if constexpr (std::is_same_v<TDirections,AllDirections>) return TStencil::Q;
        else if constexpr (std::is_same_v<TDirections,One>) return 1;
        else return TStencil::D; 
    } 
    
};

template<class TDirections=Cartesian>
struct InterfaceGradient : GradientBase<TDirections> {
    
};