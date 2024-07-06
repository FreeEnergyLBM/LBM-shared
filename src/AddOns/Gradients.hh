#pragma once
#include <iostream>

#include "../GradientStencils/GradientStencils.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "../Stencil.hh"
#include "../Template.hh"
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

    if (Geometry<Lattice>::isBulkSolid(k) || k < Lattice::HaloSize || k >= Lattice::N - Lattice::HaloSize) return;

    constexpr int numdir = TGradientStencil::template getNumberOfDirections<Stencil>();

    for (int idx = 0; idx < numdir; idx++) {
        GradientType::template get<Lattice, numdir>(k, idx) =
            mGradientStencil.template compute<TTraits, TParam>(idx, k);
    }
}

template <class TParam, class TGradientStencil = CentralXYZ>
class GradientsDirectional : public AddOnBase {
   public:
    GradientsDirectional() = default;

    GradientsDirectional(const GradientsDirectional<TGradientStencil, TParam>& other){};

    GradientsDirectional(Gradients<TGradientStencil, TParam>& other){};

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
inline void GradientsDirectional<TParam, TGradientStencil>::compute(int k) {  // Not necessary

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;
    using GradientType = typename TGradientStencil::template GradientType<TParam>;

    if (Geometry<Lattice>::isBulkSolid(k) || k < Lattice::HaloSize || k >= Lattice::N - Lattice::HaloSize) return;

    constexpr int numdir = TGradientStencil::template getNumberOfDirections<Stencil>();

    for (int idx = 0; idx < numdir; idx++) {
        for (int idx2 = 0; idx2 < numdir; idx2++) {
            GradientType::template get<Lattice, numdir, numdir>(k, idx, idx2) =
                mGradientStencil.template compute<TTraits, TParam, numdir>(idx, k, idx2);
        }
    }
}

/**
 * \brief This can calculate multiple gradients for multiple parameters
 * \tparam TParameters Tuple of parameters to calculate gradients
 * \tparam TGradientStencils Tuple of stencils to apply
 */
template <typename TParameters, typename TGradientStencils>
class GradientsMulti : public AddOnBase {
   public:
    GradientsMulti() = default;

    GradientsMulti(const GradientsMulti<TParameters, TGradientStencils>& other){};

    GradientsMulti(GradientsMulti<TParameters, TGradientStencils>& other){};

    template <class TTraits>
    inline void compute(int k);  // Perform any neccessary computations before force is computed

    inline void setWettingPrefactor(double value) {
        std::apply([value](auto&... gradient) { (gradient.setWettingPrefactor(value), ...); }, mt_ParamGradients);
    }

    inline void setInterfaceDistance(double (*distance)(int k, int idx)) {
        std::apply([distance](auto&... gradient) { (gradient.setInterfaceDistance(distance), ...); },
                   mt_ParamGradients);
    }

    inline void setInterfaceVal(double value) {
        std::apply([value](auto&... gradient) { (gradient.setInterfaceVal(value), ...); }, mt_ParamGradients);
    }

    inline void setBoundaryID(int id, bool preset = false) {
        std::apply([id, preset](auto&... gradient) { (gradient.setBoundaryID(id, preset), ...); }, mt_ParamGradients);
    }

    inline void setBoundaryID(const std::vector<int>& id, bool preset = false) {
        std::apply([id, preset](auto&... gradient) { (gradient.setBoundaryID(id, preset), ...); }, mt_ParamGradients);
    }

   private:
    using TPairs = tuple_combinations<TParameters, TGradientStencils>;
    tuple_pair_template<Gradients, TPairs> mt_ParamGradients;
};

template <typename TParameters, typename TGradientStencils>
template <class TTraits>
inline void GradientsMulti<TParameters, TGradientStencils>::compute(int k) {
    std::apply([k](auto&... gradient) { (gradient.template compute<TTraits>(k), ...); }, mt_ParamGradients);
}

/**
 * \brief This can calculate multiple gradients for multiple instances of a parameter
 * \tparam TParameter Parameter template type containing the instances to use (e.g. OrderParameter)
 * \tparam nInstances The number of instances to use
 * \tparam TGradientStencils Tuple or list of stencils to use
 */
template <template <int> typename TParameter, int nInstances, typename... TGradientStencils>
class GradientsMultiInstance
    : public GradientsMulti<int_template<TParameter, int_sequence<nInstances>>, std::tuple<TGradientStencils...>> {};

template <template <int> typename TParameter, int nInstances, typename... TGradientStencils>
class GradientsMultiInstance<TParameter, nInstances, std::tuple<TGradientStencils...>>
    : public GradientsMulti<int_template<TParameter, int_sequence<nInstances>>, std::tuple<TGradientStencils...>> {};

/**
 * \brief This can calculate a single gradient for multiple parameters
 * \tparam TGradientStencil Gradient stencil to use
 * \tparam TParameters List of Parameters to compute gradients for
 */
template <class TGradientStencil, class... TParameters>
class GradientsMultiParam : public GradientsMulti<std::tuple<TParameters...>, std::tuple<TGradientStencil>> {};

/**
 * \brief This can calculate multiple gradients for a single parameter
 * \tparam TParameter Parameter to compute gradients for
 * \tparam TGradientStencils List of gradient stencils to use
 */
template <class TParameter, class... TGradientStencils>
class GradientsMultiStencil : public GradientsMulti<std::tuple<TParameter>, std::tuple<TGradientStencils...>> {};
