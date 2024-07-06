#pragma once
#include <tuple>

#include "Collide.hh"
#include "Template.hh"

/**
 * \file Trait.hh
 * \brief Contains the DefaultTrait for LB models, and the BaseTrait object to modify it.
 */

template <class TTrait>
struct BaseTrait;

/**
 * \brief This class provides a starting point to build trait classes from.
 * \tparam TLattice The lattice the stencil will be used on.
 * \t_NumberOfComponents The number of fluid components in the simulation.
 * \details It sets a default stencil based on the number of cartesian directions and sets the collision operator to SRT
 */
template <class TLattice = void, int t_NumberOfComponents = 1>
struct DefaultTrait : BaseTrait<DefaultTrait<TLattice, t_NumberOfComponents>> {
    using Stencil =
        std::conditional_t<TLattice::NDIM == 1, D1Q3,
                           std::conditional_t<TLattice::NDIM == 2, D2Q9, D3Q19>>;  // Here, D refers to the number of
                                                                                   // cartesian dimensions

    using Boundaries = std::tuple<std::tuple<>>;

    using Processors = std::tuple<std::tuple<>>;

    using Forces = std::tuple<>;

    template <class TStencil>
    using CollisionModel = SRT<TStencil>;

    using Lattice = TLattice;

    template <class Tlattice, class TStencil>
    using DataType = DataOldNew<TLattice, TStencil>;

    static constexpr int NumberOfComponents = t_NumberOfComponents;
};

/**
 * \brief This class contains the functions for modifying a trait.
 * \tparam TTrait The starting trait that is modified.
 */
template <class TTrait>
struct BaseTrait {
    template <class... TProcessor>
    struct AddProcessor;

    template <int idx, class... TProcessor>
    struct AddProcessorIdx;

    template <class... TProcessor>
    struct AddProcessor : BaseTrait<AddProcessor<TProcessor...>> {
        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename TTrait::Boundaries;

        using Processors = std::tuple<
            tuple_cat_t<typename std::tuple_element<0, typename TTrait::Processors>::type, std::tuple<TProcessor...>>>;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <class... TProcessor>
    struct AddProcessor<std::tuple<TProcessor...>> : BaseTrait<AddProcessor<std::tuple<TProcessor...>>> {
        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename TTrait::Boundaries;

        using Processors = tuple_cat_t<typename TTrait::Processors, std::tuple<std::tuple<TProcessor...>>>;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <int idx, class... TProcessor>
    struct AddProcessorIdx : BaseTrait<AddProcessorIdx<idx, TProcessor...>> {
        static_assert(idx >= 0, "ERROR: idx must be positive");
        static_assert(idx <= std::tuple_size<typename TTrait::Processors>::value,
                      "ERROR: idx must less than or equal to the tuple size");

        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename TTrait::Boundaries;

        using Processors =
            typename add_tuple_idx<typename TTrait::Processors, idx,
                                   (std::tuple_size<typename TTrait::Processors>::value > idx), TProcessor...>::type;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <int idx, class... TProcessor>
    struct AddProcessorIdx<idx, std::tuple<TProcessor...>> : BaseTrait<AddProcessorIdx<idx, TProcessor...>> {
        static_assert(idx >= 0, "ERROR: idx must be positive");
        static_assert(idx <= std::tuple_size<typename TTrait::Processors>::value,
                      "ERROR: idx must less than or equal to the tuple size");

        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename TTrait::Boundaries;

        using Processors = typename insert_tuple_idx<typename TTrait::Processors, idx, std::tuple<TProcessor...>>::type;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <class... TProcessor>
    struct SetProcessor : BaseTrait<SetProcessor<TProcessor...>> {
        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename TTrait::Boundaries;

        using Processors = typename is_tuple<TProcessor...>::type;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <class... TForce>
    struct AddForce : BaseTrait<AddForce<TForce...>> {
        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename TTrait::Boundaries;

        using Processors = typename TTrait::Processors;

        using Forces = tuple_cat_t<typename TTrait::Forces, std::tuple<TForce...>>;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <class... TForce>
    struct SetForce : BaseTrait<SetForce<TForce...>> {
        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename TTrait::Boundaries;

        using Processors = typename TTrait::Processors;

        using Forces = std::tuple<TForce...>;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <class... TBoundary>
    struct AddBoundary;

    template <int idx, class... TBoundary>
    struct AddBoundaryIdx;

    template <class... TBoundary>
    struct AddBoundary : BaseTrait<AddBoundary<TBoundary...>> {
        using Stencil = typename TTrait::Stencil;

        using Boundaries = std::tuple<
            tuple_cat_t<typename std::tuple_element<0, typename TTrait::Boundaries>::type, std::tuple<TBoundary...>>>;

        using Processors = typename TTrait::Processors;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <class... TBoundary>
    struct AddBoundary<std::tuple<TBoundary...>> : BaseTrait<AddBoundary<std::tuple<TBoundary...>>> {
        using Stencil = typename TTrait::Stencil;

        using Boundaries = tuple_cat_t<typename TTrait::Boundaries, std::tuple<std::tuple<TBoundary...>>>;

        using Processors = typename TTrait::Processors;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <int idx, class... TBoundary>
    struct AddBoundaryIdx : BaseTrait<AddBoundaryIdx<idx, TBoundary...>> {
        static_assert(idx >= 0, "ERROR: idx must be positive");
        static_assert(idx <= std::tuple_size<typename TTrait::Boundaries>::value,
                      "ERROR: idx must less than or equal to the tuple size");

        using Stencil = typename TTrait::Stencil;

        using Boundaries =
            typename add_tuple_idx<typename TTrait::Boundaries, idx,
                                   (std::tuple_size<typename TTrait::Boundaries>::value > idx), TBoundary...>::type;

        using Processors = typename TTrait::Processors;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <int idx, class... TBoundary>
    struct AddBoundaryIdx<idx, std::tuple<TBoundary...>> : BaseTrait<AddBoundaryIdx<idx, TBoundary...>> {
        static_assert(idx >= 0, "ERROR: idx must be positive");
        static_assert(idx <= std::tuple_size<typename TTrait::Boundaries>::value,
                      "ERROR: idx must less than or equal to the tuple size");

        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename insert_tuple_idx<typename TTrait::Boundaries, idx, std::tuple<TBoundary...>>::type;

        using Processors = typename TTrait::Processors;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <class... TBoundary>
    struct SetBoundary : BaseTrait<SetBoundary<TBoundary...>> {
        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename is_tuple<TBoundary...>::type;

        using Processors = typename TTrait::Processors;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <class TStencil>
    struct SetStencil : BaseTrait<SetStencil<TStencil>> {
        using Stencil = TStencil;

        using Boundaries = typename TTrait::Boundaries;

        using Processors = typename TTrait::Processors;

        using Forces = typename TTrait::Forces;

        template <class stencil1>
        using CollisionModel = typename TTrait::template CollisionModel<stencil1>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil1>
        using DataType = typename TTrait::template DataType<TLattice, TStencil1>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <template <class> class TModel>
    struct SetCollisionOperator : BaseTrait<SetCollisionOperator<TModel>> {
        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename TTrait::Boundaries;

        using Processors = typename TTrait::Processors;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = TModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <template <class, class, bool> class TDataType, bool TSeperateStream = false>
    struct SetDataType : BaseTrait<SetDataType<TDataType>> {
        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename TTrait::Boundaries;

        using Processors = typename TTrait::Processors;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = TDataType<TLattice, TStencil, TSeperateStream>;

        static constexpr int NumberOfComponents = TTrait::NumberOfComponents;
    };

    template <int TNumberOfComponents>
    struct SetNumberOfComponents : BaseTrait<SetNumberOfComponents<TNumberOfComponents>> {
        using Stencil = typename TTrait::Stencil;

        using Boundaries = typename TTrait::Boundaries;

        using Processors = typename TTrait::Processors;

        using Forces = typename TTrait::Forces;

        template <class TStencil>
        using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

        using Lattice = typename TTrait::Lattice;

        template <class TLattice, class TStencil>
        using DataType = typename TTrait::template DataType<TLattice, TStencil>;

        static constexpr int NumberOfComponents = TNumberOfComponents;
    };
};
