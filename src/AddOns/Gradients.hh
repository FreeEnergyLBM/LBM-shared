#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "../Stencil.hh"
#include "../GradientStencils/GradientStencils.hh"

#include "AddOnBase.hh"
#include<iostream>

template<class TParam, class TGradientStencil = CentralXYZ>
class Gradients : public AddOnBase {

    public:

        Gradients() = default;

        Gradients(const Gradients<TGradientStencil,TParam>& other) {};

        Gradients(Gradients<TGradientStencil,TParam>& other) {};

        template<class TTraits>
        inline void compute(int k);

};

template<class TParam, class TGradientStencil>
template<class TTraits>
inline void Gradients<TParam, TGradientStencil>::compute(int k) { //Not necessary

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    using GradientType = typename TGradientStencil::template GradientType<TParam>;
    constexpr int numdir = TGradientStencil::template getNumberOfDirections<Stencil>();
    
    for (int component = 0 ; component < TParam::instances; component++){

        for(int idx = 0; idx < numdir; idx++) {

            GradientType::template get<Lattice,numdir>(k, component, idx) = TGradientStencil::template compute<TTraits,TParam>(idx, k, component);
            
            //std::cout<<typeid(GradientType).name()<<" "<<GradientType::template get<Lattice,numdir>(k, component, idx)<<std::endl;
        }

    }

}

template<class TGradientStencil, class ...TParam>
class GradientsMultiParam : public AddOnBase {

    public:

        GradientsMultiParam() = default;

        GradientsMultiParam(const GradientsMultiParam<TGradientStencil,TParam...>& other) {};

        GradientsMultiParam(GradientsMultiParam<TGradientStencil,TParam...>& other) {};

        template<class TTraits>
        inline void compute(int k);

    private:

        std::tuple<Gradients<TParam, TGradientStencil>...> mt_Param;

};


template<class TGradientStencil, class ...TParam>
template<class TTraits>
inline void GradientsMultiParam<TGradientStencil,TParam...>::compute(int k) { //Not necessary
    
    std::apply([k](auto&... gradient){

        (gradient.template compute<TTraits>(k), ...);

                                     }, mt_Param);

}

template<class TParam, class ...TGradientStencil>
class GradientsMultiStencil : public AddOnBase {

    public:

        GradientsMultiStencil() = default;

        GradientsMultiStencil(const GradientsMultiStencil<TParam, TGradientStencil...>& other) {};

        GradientsMultiStencil(GradientsMultiStencil<TParam, TGradientStencil...>& other) {};

        template<class TTraits>
        inline void compute(int k); //Perform any neccessary computations before force is computed

    private:

        std::tuple<Gradients<TParam, TGradientStencil>...> mt_GradientStencil;

};


template<class TParam, class ...TGradientStencil>
template<class TTraits>
inline void GradientsMultiStencil<TParam, TGradientStencil...>::compute(int k) { //Not necessary

    std::apply([k](auto&... gradient){
        (gradient.template compute<TTraits>(k), ...);
                                     }, mt_GradientStencil);

}


template<class TParam, class TGradientStencil = CentralXYZ>
class GradientsInterface : public AddOnBase {

    public:

        GradientsInterface() = default;

        GradientsInterface(const GradientsInterface<TGradientStencil,TParam>& other) {};

        GradientsInterface(GradientsInterface<TGradientStencil,TParam>& other) {};

        template<class TTraits>
        inline void compute(int k);

        inline void setInterfaceDistance(double (*distance)(int k, int idx)){

            //evalInterfaceDistance=distance;
            if constexpr (is_base_of_template<InterfaceGradient, TGradientStencil>()) mGradientStencil.setInterfaceDistance(distance);

        }

        inline void setInterfaceCondition(bool (*condition)(int k)){

            //evalInterfaceDistance=distance;
            if constexpr (is_base_of_template<InterfaceGradient, TGradientStencil>()) mGradientStencil.setInterfaceCondition(condition);

        }

        inline void setInterfaceVal(double value){

            //evalInterfaceDistance=distance;
            if constexpr (is_base_of_template<InterfaceGradient, TGradientStencil>()) mGradientStencil.setInterfaceVal(value);

        }

    private:

        TGradientStencil mGradientStencil;

        //static double defaultDistance(const double& val, int k) { return 0.5; }

        //double (*evalInterfaceDistance)(const double& val, int k) = &defaultDistance;

};

template<class TParam, class TGradientStencil>
template<class TTraits>
inline void GradientsInterface<TParam, TGradientStencil>::compute(int k) { //Not necessary

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    using GradientType = typename TGradientStencil::template GradientType<TParam>;
    constexpr int numdir = TGradientStencil::template getNumberOfDirections<Stencil>();
    
    for (int component = 0 ; component < TParam::instances; component++){

        for(int idx = 0; idx < numdir; idx++) {

            GradientType::template get<Lattice,numdir>(k, component, idx) = mGradientStencil.template compute<TTraits,TParam>(idx, k, component);
            
            //std::cout<<typeid(GradientType).name()<<" "<<GradientType::template get<Lattice,numdir>(k, component, idx)<<std::endl;
        }

    }

}