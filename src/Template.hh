#pragma once
#include <cstddef>
#include <tuple>
#include <type_traits>

/**
 * \file Template.hh
 * \brief Contains utility functions used for template metaprogramming.
 */

template <typename T>
struct remove_const_and_reference {
    using type = typename std::remove_const<typename std::remove_reference<T>::type>::type;
};

// Answer one simple question: here's a type, and a tuple. Tell me
// if the type is one of the tuples types. If so, I want it.

template <typename TWantedType, typename T>
struct is_wanted_type;

template <typename TWantedType, typename... TTypes>
struct is_wanted_type<TWantedType, std::tuple<TTypes...>> {
    static constexpr bool wanted = (std::is_same_v<TWantedType, TTypes> || ...);
};

// Ok, the ith index in the tuple, here's its std::tuple_element type.
// And TWantedElement is a tuple of all types we want to extract.
//
// Based on which way the wind blows we'll produce either a std::tuple<>
// or a std::tuple<tuple_element_t>.

template <size_t i, typename tuple_element_t, typename TWantedElement,
          bool wanted = is_wanted_type<tuple_element_t, TWantedElement>::wanted>
struct extract_type {
    template <typename TTupleType>
    inline static auto do_extract_type(TTupleType& t) {
        return std::tuple<>{};
    }
};

template <size_t i, typename tuple_element_t, typename TWantedElement>
struct extract_type<i, tuple_element_t, TWantedElement, true> {
    template <typename TTupleType>
    inline static auto do_extract_type(TTupleType& t) {
        return std::tie(std::get<i>(t));
    }
};

// And now, a simple fold expression to pull out all wanted types
// and tuple-cat them together.

template <typename TWantedElement, typename TTupleType, size_t... i>
inline auto get_type_t(TTupleType& t, std::index_sequence<i...>) {
    return std::tuple_cat(
        extract_type<i, typename std::tuple_element<i, TTupleType>::type, TWantedElement>::do_extract_type(t)...);
}

template <typename... TWantedElement, typename... types>
inline auto get_type(std::tuple<types...>& t) {
    return get_type_t<std::tuple<TWantedElement...>>(t, std::make_index_sequence<sizeof...(types)>());
}

template <class TBase, typename TTuple>
struct CheckBase;

template <class TBase, typename... TTypes>
struct CheckBase<TBase, std::tuple<TTypes...>> : std::conjunction<std::is_base_of<TBase, TTypes>...> {};

template <template <typename...> class TC, typename... Ts>
std::true_type is_base_of_template_impl(const TC<Ts...>*);

template <template <typename...> class TC>
std::false_type is_base_of_template_impl(...);

template <template <typename...> class TC, typename T>
using is_base_of_template = decltype(is_base_of_template_impl<TC>(std::declval<T*>()));

template <template <class> class TBase, typename TTuple>
struct CheckBaseTemplate;

template <template <class> class TBase, typename... TTypes>
struct CheckBaseTemplate<TBase, std::tuple<TTypes...>> : std::conjunction<is_base_of_template<TBase, TTypes>...> {};

template <typename... TInput>
using tuple_cat_t = decltype(std::tuple_cat(std::declval<TInput>()...));

template <int N, int idx, class element, class origtuple, class tuplebeforeelement>
struct tuple_replace {
    using type = typename tuple_replace<
        N - 1, idx, element, origtuple,
        tuple_cat_t<std::tuple<typename std::tuple_element<N, origtuple>::type>, tuplebeforeelement>>::type;
};

template <int N, int idx, class element, class origtuple>
struct tuple_replace<N, idx, element, origtuple, std::tuple<>> {
    using type = typename std::conditional<
        (N - 1 == idx), typename tuple_replace<N - 2, idx, element, origtuple, std::tuple<element>>::type,
        typename tuple_replace<N - 2, idx, element, origtuple,
                               std::tuple<typename std::tuple_element<N - 1, origtuple>::type>>::type>::type;
};

template <int idx, class element, class origtuple, class tuplebeforeelement>
struct tuple_replace<idx, idx, element, origtuple, tuplebeforeelement> {
    using type = typename tuple_replace<idx - 1, -999, element, origtuple,
                                        tuple_cat_t<std::tuple<element>, tuplebeforeelement>>::type;
};

template <class element, class origtuple, class tuplebeforeelement>
struct tuple_replace<-1, -999, element, origtuple, tuplebeforeelement> {
    using type = tuplebeforeelement;
};

template <class element, class origtuple, class tuplebeforeelement>
struct tuple_replace<-1, 0, element, origtuple, tuplebeforeelement> {
    using type = tuplebeforeelement;
};

template <class element, class origtuple, class tuplebeforeelement>
struct tuple_replace<0, -999, element, origtuple, tuplebeforeelement> {
    using type = tuple_cat_t<std::tuple<typename std::tuple_element<0, origtuple>::type>, tuplebeforeelement>;
};

template <int N, int idx, class element, class origtuple, class tuplebeforeelement>
struct tuple_insert {
    using type = typename tuple_insert<
        N - 1, idx, element, origtuple,
        tuple_cat_t<std::tuple<typename std::tuple_element<N, origtuple>::type>, tuplebeforeelement>>::type;
};

template <int N, int idx, class element, class origtuple>
struct tuple_insert<N, idx, element, origtuple, std::tuple<>> {
    using type = typename std::conditional<
        (idx >= N),
        typename tuple_insert<
            N - 2, idx, element, origtuple,
            tuple_cat_t<std::tuple<typename std::tuple_element<N - 1, origtuple>::type>, std::tuple<element>>>::type,
        typename std::conditional<
            (N - 1 == idx),
            typename tuple_insert<N - 2, idx, element, origtuple,
                                  tuple_cat_t<std::tuple<element>,
                                              std::tuple<typename std::tuple_element<N - 1, origtuple>::type>>>::type,
            typename tuple_insert<N - 2, idx, element, origtuple,
                                  std::tuple<typename std::tuple_element<N - 1, origtuple>::type>>::type>::type>::type;
};

template <int idx, class element, class origtuple, class tuplebeforeelement>
struct tuple_insert<idx, idx, element, origtuple, tuplebeforeelement> {
    using type = typename tuple_insert<
        idx - 1, -999, element, origtuple,
        tuple_cat_t<tuple_cat_t<std::tuple<element>, std::tuple<typename std::tuple_element<idx, origtuple>::type>>,
                    tuplebeforeelement>>::type;
};

template <class element, class origtuple, class tuplebeforeelement>
struct tuple_insert<0, -999, element, origtuple, tuplebeforeelement> {
    using type = tuple_cat_t<std::tuple<typename std::tuple_element<0, origtuple>::type>, tuplebeforeelement>;
};

template <class element, class origtuple, class tuplebeforeelement>
struct tuple_insert<-1, -999, element, origtuple, tuplebeforeelement> {
    using type = tuplebeforeelement;
};

template <class element, class origtuple, class tuplebeforeelement>
struct tuple_insert<-1, 0, element, origtuple, tuplebeforeelement> {
    using type = tuplebeforeelement;
};

template <typename... T>
struct is_tuple {
    using type = std::tuple<std::tuple<T...>>;
};

template <typename... T1, typename... T2>
struct is_tuple<std::tuple<T1...>, T2...> {
    using type = std::tuple<std::tuple<T1...>, T2...>;
};

template <typename tups, int idx, bool intuple, typename... T>
struct add_tuple_idx;

template <typename tups, int idx, typename... T>
struct insert_tuple_idx;

template <typename tups, int idx, typename... T>
struct add_tuple_idx<tups, idx, true, T...> {
    using type = typename tuple_replace<std::tuple_size<tups>::value, idx,
                                        tuple_cat_t<typename std::tuple_element<idx, tups>::type, std::tuple<T...>>,
                                        tups, std::tuple<>>::type;
};

template <typename tups, int idx, typename... T>
struct add_tuple_idx<tups, idx, false, T...> {
    using type = tuple_cat_t<tups, std::tuple<std::tuple<T...>>>;
};

template <typename tups, int idx, typename... T>
struct insert_tuple_idx<tups, idx, std::tuple<T...>> {
    using type = typename tuple_insert<std::tuple_size<tups>::value, idx, std::tuple<T...>, tups, std::tuple<>>::type;
};

template <std::size_t i, class TTuple, std::size_t... is>
constexpr auto element_as_tuple(const TTuple tuple, std::index_sequence<is...>) {
    if constexpr (!(std::is_same_v<std::tuple_element_t<i, TTuple>, std::tuple_element_t<is, TTuple>> || ...))
        return std::make_tuple(std::get<i>(tuple));
    else
        return std::make_tuple();
}

template <class TTuple, std::size_t... is>
constexpr auto make_tuple_unique(const TTuple tuple, std::index_sequence<is...>) {
    return std::tuple_cat(element_as_tuple<is>(tuple, std::make_index_sequence<is>{})...);
}

template <class... Tuples>
constexpr auto make_tuple_unique(const Tuples... tuples) {
    auto all = std::tuple_cat(tuples...);
    constexpr auto size = std::tuple_size_v<decltype(all)>;
    return make_tuple_unique(all, std::make_index_sequence<size>{});
}

template <typename T, typename TTuple>
struct has_type;

template <typename T>
struct has_type<T, std::tuple<>> : std::false_type {};

template <typename T, typename U, typename... Ts>
struct has_type<T, std::tuple<U, Ts...>> : has_type<T, std::tuple<Ts...>> {};

template <typename T, typename... Ts>
struct has_type<T, std::tuple<T, Ts...>> : std::true_type {};

struct Cartesian {};
struct AllDirections {};
struct One {};

template <typename TKey, typename TValue>
struct kv {
    using Key = TKey;

    static TValue Value;
};

template <typename TKey, typename TValue>
TValue kv<TKey, TValue>::Value;

template <typename...>
struct ct_map;

template <>
struct ct_map<> {
    template <typename>
    struct keyexists {
        static constexpr bool exists = false;
    };

    template <typename T>
    struct get {
        static inline constexpr int noKey() {
            if constexpr (!sizeof(T)) {
                static_assert(!!sizeof(T), "Key does not exist in map.");
            }
            return 0;
        }

        static constexpr auto val = noKey();
    };
};

template <typename TKey, typename TValue, typename... TRest>
struct ct_map<kv<TKey, TValue>, TRest...> {
    template <typename TKKey>
    struct keyexists {
        static constexpr bool exists = (std::is_same<typename std::remove_reference<TKKey>::type,
                                                     typename std::remove_reference<TKey>::type>::value)
                                           ? true
                                           : ct_map<TRest...>::template keyexists<TKKey>::exists;
    };

    template <typename TKKey>
    struct get {
        static inline constexpr auto& findVal() {
            if constexpr (sizeof...(TRest) != 0) {
                if constexpr (std::is_same<TKKey, TKey>::value) {
                    return kv<TKey, TValue>::Value;
                } else {
                    static_assert(
                        sizeof...(TRest) != 0 || std::is_same<typename std::remove_reference<TKKey>::type,
                                                              typename std::remove_reference<TKey>::type>::value,
                        "Key does not exist in map.");
                    return ct_map<TRest...>::template get<TKKey>::val;
                }

            } else {
                return kv<TKey, TValue>::Value;
            }
        }

        static constexpr auto& val = findVal();
    };
};

template <typename TKey, typename TValue>
struct kv_types {
    using Key = TKey;

    using Type = TValue;
};

template <typename...>
struct ct_map_types;

template <>
struct ct_map_types<> {
    template <typename>
    struct keyexists {
        static constexpr bool exists = false;
    };

    template <typename T>
    struct get {
        static inline constexpr int noKey() {
            if constexpr (!sizeof(T)) {
                static_assert(!!sizeof(T), "Key does not exist in map.");
            }
            return 0;
        }

        using TValue = decltype(noKey());

        TValue val = noKey();
    };
};

template <typename TKey, typename TValue, typename... TRest>
struct ct_map_types<kv<TKey, TValue>, TRest...> {
    template <typename TKKey>
    struct keyexists {
        static constexpr bool exists = (std::is_same<typename std::remove_reference<TKKey>::type,
                                                     typename std::remove_reference<TKey>::type>::value)
                                           ? true
                                           : ct_map_types<TRest...>::template keyexists<TKKey>::exists;
    };

    template <typename TKKey>
    struct get {
        static inline constexpr auto findVal() {
            static_assert(sizeof...(TRest) != 0 || std::is_same<typename std::remove_reference<TKKey>::type,
                                                                typename std::remove_reference<TKey>::type>::value,
                          "Key does not exist in map.");
            if constexpr (sizeof...(TRest) != 0) {
                if constexpr (std::is_same<TKKey, TKey>::value) {
                    typename kv_types<TKey, TValue>::Type val = {};
                    return val;
                } else {
                    typename ct_map_types<TRest...>::template get<TKKey>::Type val = {};
                    return val;
                }
            } else {
                typename kv_types<TKey, TValue>::Type val = {};
                return val;
            }
        }

        using Type = decltype(findVal());

        Type val = findVal();
    };
};

template <typename T, typename C, int I>
struct tuple_index_r;

template <typename H, typename... R, typename C, int I>
struct tuple_index_r<std::tuple<H, R...>, C, I>
    : public std::conditional<std::is_same<C, typename std::remove_reference<H>::type>::value,
                              std::integral_constant<int, I>, tuple_index_r<std::tuple<R...>, C, I + 1>>::type {};

template <typename C, int I>
struct tuple_index_r<std::tuple<>, C, I> : std::integral_constant<int, -1> {};

template <typename T, typename C>
struct tuple_index_of : public std::integral_constant<int, tuple_index_r<T, C, 0>::value> {};

template <typename T, typename C, int I>
struct tuple_tuple_index_r;

template <typename H, typename... R, typename C, int I>
struct tuple_tuple_index_r<std::tuple<H, R...>, C, I> {
    static constexpr int outeridx = tuple_index_of<H, C>::value;
    using idx =
        typename std::conditional<(outeridx >= 0),
                                  std::tuple<std::integral_constant<int, I>, std::integral_constant<int, outeridx>>,
                                  typename tuple_tuple_index_r<std::tuple<R...>, C, I + 1>::idx>::type;
};

template <typename C, int I>
struct tuple_tuple_index_r<std::tuple<>, C, I> {
    using idx = int;
};

template <typename T, typename C>
struct tuple_tuple_index_of : public tuple_tuple_index_r<T, typename std::remove_reference<C>::type, 0> {};

// Check if a type is a std::tuple
template <typename T>
struct is_std_tuple : std::false_type {};

template <typename... Ts>
struct is_std_tuple<std::tuple<Ts...>> : std::true_type {};

// Get all combinations of two tuples as a tuple of pairs
template <typename Tuple1, typename Tuple2>
struct tuple_combinations_impl;

template <typename... T1s, typename... T2s>
struct tuple_combinations_impl<std::tuple<T1s...>, std::tuple<T2s...>> {
    // Get combination of one T1 with all T2s
    template <typename T1>
    struct combine {
        using type = std::tuple<std::pair<T1, T2s>...>;
    };

    using type = tuple_cat_t<typename combine<T1s>::type...>;
};

template <typename Tuple1, typename Tuple2>
using tuple_combinations = typename tuple_combinations_impl<Tuple1, Tuple2>::type;

// Metafunction to assign the templates of a type using an std::pair
template <template <typename, typename> typename TargetType, typename Pair>
struct pair_template_impl;

template <template <typename, typename> typename TargetType, typename A, typename B>
struct pair_template_impl<TargetType, std::pair<A, B>> {
    using type = TargetType<A, B>;
};

template <template <typename, typename> typename TargetType, typename Pair>
using pair_template = typename pair_template_impl<TargetType, Pair>::type;

// Metafunction to apply pair_template across a tuple of pairs
template <template <typename, typename> typename TargetType, typename Tuple>
struct tuple_pair_template_impl;

template <template <typename, typename> typename TargetType, typename... Pairs>
struct tuple_pair_template_impl<TargetType, std::tuple<Pairs...>> {
    using type = std::tuple<pair_template<TargetType, Pairs>...>;
};

template <template <typename, typename> typename TargetType, typename Tuple>
using tuple_pair_template = typename tuple_pair_template_impl<TargetType, Tuple>::type;

// Metafunction to create a tuple of types using a tuple of ints
template <template <int> typename TargetType, typename Integers>
struct int_template_impl;

template <template <int> typename TargetType, typename T, T... Is>
struct int_template_impl<TargetType, std::integer_sequence<T, Is...>> {
    using type = std::tuple<TargetType<Is>...>;
};

template <template <int> typename TargetType, typename Integers>
using int_template = typename int_template_impl<TargetType, Integers>::type;

// Generate an integer_sequence up to N
template <int N, int... Is>
struct int_sequence_impl : int_sequence_impl<N - 1, N - 1, Is...> {};

template <int... Is>
struct int_sequence_impl<0, Is...> {
    using type = std::integer_sequence<int, Is...>;
};

template <int N, int... Is>
using int_sequence = typename int_sequence_impl<N, Is...>::type;
