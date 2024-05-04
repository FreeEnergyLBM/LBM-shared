#pragma once
#include <iostream>

#include "../GradientStencils/GradientStencils.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "../Stencil.hh"
#include "AddOnBase.hh"

template <class TParam, class TGradientStencil = CentralXYZ>
class Gradients : public AddOnBase {
   public:
    Gradients() = default;

    Gradients(const Gradients<TGradientStencil, TParam>& other){};

    Gradients(Gradients<TGradientStencil, TParam>& other){};

    template <class TTraits>
    inline void compute(int k);

    inline void setInterfaceDistance(double (*distance)(int k, int idx)) {
        mGradientStencil.setInterfaceDistance(distance);
    }

    inline void setInterfaceVal(double value) { mGradientStencil.setInterfaceVal(value); }

    inline void setWettingPrefactor(double value) { mGradientStencil.setPrefactor(value); }

    inline void setBoundaryID(int id, bool preset = false) { mGradientStencil.setBoundaryID(id, preset); }

    inline void setBoundaryID(const std::vector<int>& id, bool preset = false) {
        mGradientStencil.setBoundaryID(id, preset);
    }

   private:
    TGradientStencil mGradientStencil;
};

template <class TParam, class TGradientStencil>
template <class TTraits>
inline void Gradients<TParam, TGradientStencil>::compute(int k) {  // Not necessary

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    using GradientType = typename TGradientStencil::template GradientType<TParam>;

    if (Geometry<Lattice>::isBulkSolid(k)) return;

    constexpr int numdir = TGradientStencil::template getNumberOfDirections<Stencil>();

    for (int component = 0; component < TParam::instances; component++) {
        for (int idx = 0; idx < numdir; idx++) {
            GradientType::template get<Lattice, numdir>(k, component, idx) =
                mGradientStencil.template compute<TTraits, TParam>(idx, k, component);
        }
    }
}

template <class TGradientStencil, class... TParam>
class GradientsMultiParam : public AddOnBase {
   public:
    GradientsMultiParam() = default;

    GradientsMultiParam(const GradientsMultiParam<TGradientStencil, TParam...>& other){};

    GradientsMultiParam(GradientsMultiParam<TGradientStencil, TParam...>& other){};

    template <class TTraits>
    inline void compute(int k);

    inline void setWettingPrefactor(double value) {
        std::apply([value](auto&... gradient) { (gradient.setWettingPrefactor(value), ...); }, mt_Param);
    }

    inline void setInterfaceDistance(double (*distance)(int k, int idx)) {
        std::apply([distance](auto&... gradient) { (gradient.setInterfaceDistance(distance), ...); }, mt_Param);
    }

    inline void setInterfaceVal(double value) {
        std::apply([value](auto&... gradient) { (gradient.setInterfaceVal(value), ...); }, mt_Param);
    }

    inline void setBoundaryID(int id, bool preset = false) {
        std::apply([id, preset](auto&... gradient) { (gradient.setBoundaryID(id, preset), ...); }, mt_Param);
    }

    inline void setBoundaryID(const std::vector<int>& id, bool preset = false) {
        std::apply([id, preset](auto&... gradient) { (gradient.setBoundaryID(id, preset), ...); }, mt_Param);
    }

   private:
    std::tuple<Gradients<TParam, TGradientStencil>...> mt_Param;
};

template <class TGradientStencil, class... TParam>
template <class TTraits>
inline void GradientsMultiParam<TGradientStencil, TParam...>::compute(int k) {  // Not necessary

    std::apply([k](auto&... gradient) { (gradient.template compute<TTraits>(k), ...); }, mt_Param);
}

template <class TParam, class... TGradientStencil>
class GradientsMultiStencil : public AddOnBase {
   public:
    GradientsMultiStencil() = default;

    GradientsMultiStencil(const GradientsMultiStencil<TParam, TGradientStencil...>& other){};

    GradientsMultiStencil(GradientsMultiStencil<TParam, TGradientStencil...>& other){};

    template <class TTraits>
    inline void compute(int k);  // Perform any neccessary computations before force is computed

    inline void setWettingPrefactor(double value) {
        std::apply([value](auto&... gradient) { (gradient.setWettingPrefactor(value), ...); }, mt_GradientStencil);
    }

    inline void setInterfaceDistance(double (*distance)(int k, int idx)) {
        std::apply([distance](auto&... gradient) { (gradient.setInterfaceDistance(distance), ...); },
                   mt_GradientStencil);
    }

    inline void setInterfaceVal(double value) {
        std::apply([value](auto&... gradient) { (gradient.setInterfaceVal(value), ...); }, mt_GradientStencil);
    }

    inline void setBoundaryID(int id, bool preset = false) {
        std::apply([id, preset](auto&... gradient) { (gradient.setBoundaryID(id, preset), ...); }, mt_GradientStencil);
    }

    inline void setBoundaryID(const std::vector<int>& id, bool preset = false) {
        std::apply([id, preset](auto&... gradient) { (gradient.setBoundaryID(id, preset), ...); }, mt_GradientStencil);
    }

   private:
    std::tuple<Gradients<TParam, TGradientStencil>...> mt_GradientStencil;
};

template <class TParam, class... TGradientStencil>
template <class TTraits>
inline void GradientsMultiStencil<TParam, TGradientStencil...>::compute(int k) {  // Not necessary

    std::apply([k](auto&... gradient) { (gradient.template compute<TTraits>(k), ...); }, mt_GradientStencil);
}