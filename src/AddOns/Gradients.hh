#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "../Stencil.hh"
#include "../GradientStencils/GradientStencils.hh"

#include "AddOnBase.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<class param, class gradientstencil=CentralXYZ>
class Gradients : public AddOnBase {

    public:

        Gradients() = default;

        Gradients(const Gradients<gradientstencil,param>& other) {};

        Gradients(Gradients<gradientstencil,param>& other) {};

        template<class traits>
        inline void compute(const int k); //Perform any neccessary computations before force is computed

        template<int>
        void test(){}

};

template<class param, class gradientstencil>
template<class traits>
inline void Gradients<param, gradientstencil>::compute(const int k) { //Not necessary

    for (int component = 0 ; component<param::instances; component++){
        for(int idx = 0; idx <gradientstencil::template getNumberOfDirections<typename traits::Stencil>(); idx++) {
            //test<gradientstencil::template GradientType<param>::template get<typename traits::Lattice,gradientstencil::template getNumberOfDirections<typename traits::Stencil>()>>();
            gradientstencil::template GradientType<param>::template get<typename traits::Lattice,gradientstencil::template getNumberOfDirections<typename traits::Stencil>()>(k,component,idx) = gradientstencil::template compute<traits,param>(idx, k, component);
            
        }
    }

}

template<class gradientstencil, class ...param>
class GradientsMultiParam : public AddOnBase {

    public:

        GradientsMultiParam() = default;

        GradientsMultiParam(const GradientsMultiParam<gradientstencil,param...>& other) {};

        GradientsMultiParam(GradientsMultiParam<gradientstencil,param...>& other) {};

        template<class traits>
        inline void compute(const int k); //Perform any neccessary computations before force is computed

    private:

        std::tuple<Gradients<param,gradientstencil>...> mt_Param;

};


template<class gradientstencil, class ...param>
template<class traits>
inline void GradientsMultiParam<gradientstencil,param...>::compute(const int k) { //Not necessary
    
    std::apply([k](auto&... gradient){
        (gradient.template compute<traits>(k),...);
    },mt_Param);

}

template<class param, class ...gradientstencil>
class GradientsMultiStencil : public AddOnBase {

    public:

        GradientsMultiStencil() = default;

        GradientsMultiStencil(const GradientsMultiStencil<param,gradientstencil...>& other) {};

        GradientsMultiStencil(GradientsMultiStencil<param,gradientstencil...>& other) {};

        template<class traits>
        inline void compute(const int k); //Perform any neccessary computations before force is computed

    private:

        std::tuple<Gradients<param,gradientstencil>...> mt_GradientStencil;

};


template<class param, class ...gradientstencil>
template<class traits>
inline void GradientsMultiStencil<param, gradientstencil...>::compute(const int k) { //Not necessary

    std::apply([k](auto&... gradient){
        (gradient.template compute<traits>(k),...);
    },mt_GradientStencil);

}