#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "../Stencil.hh"
#include "../GradientStencils/GradientStencils.hh"

#include "AddOnBase.hh"
#include<iostream>

template<class T_param, class T_gradientstencil = CentralXYZ>
class Gradients : public AddOnBase {

    public:

        Gradients() = default;

        Gradients(const Gradients<T_gradientstencil,T_param>& other) {};

        Gradients(Gradients<T_gradientstencil,T_param>& other) {};

        template<class T_traits>
        inline void compute(int k);

};

template<class T_param, class T_gradientstencil>
template<class T_traits>
inline void Gradients<T_param, T_gradientstencil>::compute(int k) { //Not necessary

    using Lattice = typename T_traits::Lattice;
    using Stencil = typename T_traits::Stencil;
    using GradientType = typename T_gradientstencil::template GradientType<T_param>;
    constexpr int numdir = T_gradientstencil::template getNumberOfDirections<Stencil>();

    for (int component = 0 ; component < T_param::instances; component++){

        for(int idx = 0; idx < numdir; idx++) {

            GradientType::template get<Lattice,numdir>(k, component, idx) = T_gradientstencil::template compute<T_traits,T_param>(idx, k, component);
            
        }

    }

}

template<class T_gradientstencil, class ...T_param>
class GradientsMultiParam : public AddOnBase {

    public:

        GradientsMultiParam() = default;

        GradientsMultiParam(const GradientsMultiParam<T_gradientstencil,T_param...>& other) {};

        GradientsMultiParam(GradientsMultiParam<T_gradientstencil,T_param...>& other) {};

        template<class T_traits>
        inline void compute(int k);

    private:

        std::tuple<Gradients<T_param, T_gradientstencil>...> mt_Param;

};


template<class T_gradientstencil, class ...T_param>
template<class T_traits>
inline void GradientsMultiParam<T_gradientstencil,T_param...>::compute(int k) { //Not necessary
    
    std::apply([k](auto&... gradient){

        (gradient.template compute<T_traits>(k), ...);

                                     }, mt_Param);

}

template<class T_param, class ...T_gradientstencil>
class GradientsMultiStencil : public AddOnBase {

    public:

        GradientsMultiStencil() = default;

        GradientsMultiStencil(const GradientsMultiStencil<T_param, T_gradientstencil...>& other) {};

        GradientsMultiStencil(GradientsMultiStencil<T_param, T_gradientstencil...>& other) {};

        template<class T_traits>
        inline void compute(int k); //Perform any neccessary computations before force is computed

    private:

        std::tuple<Gradients<T_param, T_gradientstencil>...> mt_GradientStencil;

};


template<class T_param, class ...T_gradientstencil>
template<class T_traits>
inline void GradientsMultiStencil<T_param, T_gradientstencil...>::compute(int k) { //Not necessary

    std::apply([k](auto&... gradient){
        (gradient.template compute<T_traits>(k), ...);
                                     }, mt_GradientStencil);

}