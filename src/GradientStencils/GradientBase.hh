#pragma once
#include "../Service.hh"

template<class Directions=Cartesian>
struct GradientBase{

    template<class stencil>
    inline static constexpr int getNumberOfDirections(){
        if constexpr (std::is_same_v<Directions,Cartesian>) return stencil::D; 
        else if constexpr (std::is_same_v<Directions,AllDirections>) return stencil::Q;
        else if constexpr (std::is_same_v<Directions,One>) return 1;
        else return stencil::D; 
    } 
    
};