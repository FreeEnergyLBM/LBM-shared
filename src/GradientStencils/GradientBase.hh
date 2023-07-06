#pragma once
#include "../Service.hh"

template<class T_Directions=Cartesian>
struct GradientBase{

    template<class T_stencil>
    inline static constexpr int getNumberOfDirections(){
        if constexpr (std::is_same_v<T_Directions,Cartesian>) return T_stencil::D; 
        else if constexpr (std::is_same_v<T_Directions,AllDirections>) return T_stencil::Q;
        else if constexpr (std::is_same_v<T_Directions,One>) return 1;
        else return T_stencil::D; 
    } 
    
};