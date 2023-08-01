#pragma once
#include <map>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <any>
//#include<charconv>
#include <memory>
#include <cassert>
#include <complex>
#include <cstdint>
#include <type_traits>
#include <cstddef>
#include <utility>
#include <tuple>
#include <type_traits>
#include <iostream>
#ifdef MPIPARALLEL
#include <mpi.h>
#endif
#include "Stencil.hh"
#include "Mpi.hh"

//Service.hh: This will contain some commonly used functions with various uses.


template<class TLattice>
inline int computeXGlobal(const int k) //Compute X direction from a given k, the convention in this code is that
                          //k will iterate over the z direction first, then increment y by 1 once it reaches LZ,
                          //then repeat the iteration over z. Once it reaches LY x will be incremented and this
                          //process continues
{  

  return TLattice::LXMPIOffset+int((k-TLattice::HaloSize)/(float) (TLattice::LZ*TLattice::LY));


}

inline int computeX(const int& LY,const int& LZ,const int k) //Compute X direction from a given k, the convention in this code is that
                          //k will iterate over the z direction first, then increment y by 1 once it reaches LZ,
                          //then repeat the iteration over z. Once it reaches LY x will be incremented and this
                          //process continues
{  

  return int(k/(float) (LZ*LY));


}

inline int computeY(const int& LY,const int& LZ,const int k) //Compute Y direction from a given k, this uses information from the X direction
{  

  return int((k-computeX(LY,LZ,k)*LZ*LY)/(float) LZ);

}

inline int computeZ(const int& LY,const int& LZ,const int k) //Compute Y direction from a given k, this uses information from the X and Y directions
{  

  return k-computeX(LY,LZ,k)*LZ*LY-computeY(LY,LZ,k)*LZ;

}

template<class TForce>
typename TForce::Method getMethod(TForce& f){
    return std::declval<typename TForce::Method>();
}

template<class TForce>
typename TForce::Prefactor getForcePrefactor(TForce& f){
    return std::declval<typename TForce::Prefactor>();
}

void print() {
    if (mpi.rank != 0) return;
    std::cout << std::endl;
}

template <typename T, typename ... TArgs>
void print(std::vector<T> first, TArgs ... args) {
    if (mpi.rank != 0) return;
    for (auto elem : first) {
      std::cout << elem << " ";
    }
    print(args...);
}

template <typename T, typename ... TArgs>
void print(T first, TArgs ... args) {
    if (mpi.rank != 0) return;
    std::cout << first << " ";
    print(args...);
}

void printAll() {
    std::cout << std::endl;
}

template <typename T, typename ... TArgs>
void printAll(std::vector<T> first, TArgs ... args) {
    for (auto elem : first) {
      std::cout << elem << " ";
    }
    printAll(args...);
}

template <typename T, typename ... TArgs>
void printAll(T first, TArgs ... args) {
    std::cout << first << " ";
    printAll(args...);
}

template<typename T>
struct remove_const_and_reference{
  using type = typename std::remove_const<typename std::remove_reference<T>::type>::type;
};

/**\fn      mpi_get_type
 * \brief   Small template function to return the correct MPI_DATATYPE
 *          data type need for an MPI message as a constexpr at compile time
 *          https://www.mpich.org/static/docs/latest/www3/Constants.html
 *          Call in a template function with mpi_get_type<T>()
 * 
 * \tparam  T   The C++ data type used in the MPI function
 * \return  The MPI_Datatype belonging to the template C++ data type T
*/
#ifdef MPIPARALLEL
template <typename T>
[[nodiscard]] constexpr MPI_Datatype mpi_get_type() noexcept {

  MPI_Datatype mpi_type = MPI_DATATYPE_NULL;
    
  if constexpr (std::is_same_v<T, char>) {
    mpi_type = MPI_CHAR;
  }
  else if constexpr (std::is_same_v<T, signed char>) {
    mpi_type = MPI_SIGNED_CHAR;
  }
  else if constexpr (std::is_same_v<T, unsigned char>) {
    mpi_type = MPI_UNSIGNED_CHAR;
  }
  else if constexpr (std::is_same_v<T, wchar_t>) {
    mpi_type = MPI_WCHAR;
  }
  else if constexpr (std::is_same_v<T, signed short>) {
    mpi_type = MPI_SHORT;
  }
  else if constexpr (std::is_same_v<T, unsigned short>) {
    mpi_type = MPI_UNSIGNED_SHORT;
  }
  else if constexpr (std::is_same_v<T, signed int>) {
    mpi_type = MPI_INT;
  }
  else if constexpr (std::is_same_v<T, unsigned int>) {
    mpi_type = MPI_UNSIGNED;
  }
  else if constexpr (std::is_same_v<T, signed long int>) {
     mpi_type = MPI_LONG;
  }
  else if constexpr (std::is_same_v<T, unsigned long int>) {
    mpi_type = MPI_UNSIGNED_LONG;
  }
  else if constexpr (std::is_same_v<T, signed long long int>) {
    mpi_type = MPI_LONG_LONG;
  }
  else if constexpr (std::is_same_v<T, unsigned long long int>) {
    mpi_type = MPI_UNSIGNED_LONG_LONG;
  }
  else if constexpr (std::is_same_v<T, float>) {
    mpi_type = MPI_FLOAT;
  }
  else if constexpr (std::is_same_v<T, double>) {
    mpi_type = MPI_DOUBLE;
  }
  else if constexpr (std::is_same_v<T, long double>) {
    mpi_type = MPI_LONG_DOUBLE;
  }
  else if constexpr (std::is_same_v<T, std::int8_t>) {
    mpi_type = MPI_INT8_T;
  }
  else if constexpr (std::is_same_v<T, std::int16_t>) {
    mpi_type = MPI_INT16_T;
  }
  else if constexpr (std::is_same_v<T, std::int32_t>) {
    mpi_type = MPI_INT32_T;
  }
  else if constexpr (std::is_same_v<T, std::int64_t>) {
    mpi_type = MPI_INT64_T;
  }
  else if constexpr (std::is_same_v<T, std::uint8_t>) {
    mpi_type = MPI_UINT8_T;
  }
  else if constexpr (std::is_same_v<T, std::uint16_t>) {
    mpi_type = MPI_UINT16_T;
  }
  else if constexpr (std::is_same_v<T, std::uint32_t>) {
    mpi_type = MPI_UINT32_T;
  }
  else if constexpr (std::is_same_v<T, std::uint64_t>) {
    mpi_type = MPI_UINT64_T;
  }
  else if constexpr (std::is_same_v<T, bool>) {
    mpi_type = MPI_C_BOOL;
  }
  else if constexpr (std::is_same_v<T, std::complex<float>>) {
    mpi_type = MPI_C_COMPLEX;
  }
  else if constexpr (std::is_same_v<T, std::complex<double>>) {
    mpi_type = MPI_C_DOUBLE_COMPLEX;
  }
  else if constexpr (std::is_same_v<T, std::complex<long double>>) {
    mpi_type = MPI_C_LONG_DOUBLE_COMPLEX;
  }
	
  assert(mpi_type != MPI_DATATYPE_NULL);

  return mpi_type;

}
#endif


// Answer one simple question: here's a type, and a tuple. Tell me
// if the type is one of the tuples types. If so, I want it.

template<typename TWantedType, typename T> struct is_wanted_type;

template<typename TWantedType, typename ...TTypes>
struct is_wanted_type<TWantedType, std::tuple<TTypes...>> {

    static constexpr bool wanted = (std::is_same_v<TWantedType, TTypes> || ...);

};

// Ok, the ith index in the tuple, here's its std::tuple_element type.
// And TWantedElement is a tuple of all types we want to extract.
//
// Based on which way the wind blows we'll produce either a std::tuple<>
// or a std::tuple<tuple_element_t>.

template<size_t i, typename tuple_element_t,
         typename TWantedElement,
         bool wanted = is_wanted_type<tuple_element_t, TWantedElement>::wanted>
struct extract_type {

    template<typename TTupleType>
    inline static auto do_extract_type(TTupleType &t) {

        return std::tuple<>{};

    }

};


template<size_t i, typename tuple_element_t, typename TWantedElement>
struct extract_type<i, tuple_element_t, TWantedElement, true> {

    template<typename TTupleType>
    inline static auto do_extract_type(TTupleType &t) {

        return std::tie(std::get<i>(t));

    }

};

// And now, a simple fold expression to pull out all wanted types
// and tuple-cat them together.

template<typename TWantedElement, typename TTupleType, size_t ...i>
inline auto get_type_t(TTupleType &t, std::index_sequence<i...>) {

    return std::tuple_cat(extract_type<i,
                          typename std::tuple_element<i, TTupleType>::type,
                          TWantedElement>::do_extract_type(t)...);

}


template<typename ...TWantedElement, typename ...types>
inline auto get_type(std::tuple<types...> &t) {

    return get_type_t<std::tuple<TWantedElement...>>(t, std::make_index_sequence<sizeof...(types)>());
        
}

template<class TBase, typename TTuple>
struct CheckBase;

template<class TBase, typename... TTypes>
struct CheckBase<TBase, std::tuple<TTypes...>> : std::conjunction<std::is_base_of<TBase,TTypes>...> {};


template <template <typename...> class TC, typename...Ts>
std::true_type is_base_of_template_impl(const TC<Ts...>*);

template <template <typename...> class TC>
std::false_type is_base_of_template_impl(...);

template <template <typename...> class TC,typename T>
using is_base_of_template = decltype(is_base_of_template_impl<TC>(std::declval<T*>()));

template<template<class> class TBase, typename TTuple>
struct CheckBaseTemplate;

template<template<class> class TBase, typename... TTypes>
struct CheckBaseTemplate<TBase, std::tuple<TTypes...>> : std::conjunction<is_base_of_template<TBase,TTypes>...> {};

template<typename ... TInput>
using tuple_cat_t=
decltype(std::tuple_cat(
    std::declval<TInput>()...
));

template<class TTrait>
struct BaseTrait{

  template<class... TPreProcessor>
  struct AddPreProcessor : BaseTrait<AddPreProcessor<TPreProcessor...>>  {

    using Stencil = typename TTrait::Stencil;

    using Boundaries = typename TTrait::Boundaries;

    using PreProcessors = tuple_cat_t<typename TTrait::PreProcessors, std::tuple<TPreProcessor...>>;

    using PostProcessors = typename TTrait::PostProcessors;

    using Forces = typename TTrait::Forces;

    template<class TStencil>
    using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

    using Lattice = typename TTrait::Lattice;

    static constexpr int NumberOfComponents = TTrait::NumberOfComponents;

  };

  template<class... TPreProcessor>
  struct SetPreProcessor : BaseTrait<SetPreProcessor<TPreProcessor...>> {

    using Stencil = typename TTrait::Stencil;

    using Boundaries = typename TTrait::Boundaries;

    using PreProcessors = std::tuple<TPreProcessor...>;

    using PostProcessors = typename TTrait::PostProcessors;

    using Forces = typename TTrait::Forces;

    template<class TStencil>
    using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

    using Lattice = typename TTrait::Lattice;

    static constexpr int NumberOfComponents = TTrait::NumberOfComponents;

  };

  template<class... TPostProcessor>
  struct AddPostProcessor : BaseTrait<AddPostProcessor<TPostProcessor...>>  {

    using Stencil = typename TTrait::Stencil;

    using Boundaries = typename TTrait::Boundaries;

    using PreProcessors = typename TTrait::PreProcessors;

    using PostProcessors = tuple_cat_t<typename TTrait::PostProcessors, std::tuple<TPostProcessor...>>;

    using Forces = typename TTrait::Forces;

    template<class TStencil>
    using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

    using Lattice = typename TTrait::Lattice;

    static constexpr int NumberOfComponents = TTrait::NumberOfComponents;

  };

  template<class... TPostProcessor>
  struct SetPostProcessor : BaseTrait<SetPostProcessor<TPostProcessor...>> {

    using Stencil = typename TTrait::Stencil;

    using Boundaries = typename TTrait::Boundaries;

    using PreProcessors = typename TTrait::PreProcessors;

    using PostProcessors = std::tuple<TPostProcessor...>;

    using Forces = typename TTrait::Forces;

    template<class TStencil>
    using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

    using Lattice = typename TTrait::Lattice;

    static constexpr int NumberOfComponents = TTrait::NumberOfComponents;

  };
  
  template<class... TForce>
  struct AddForce : BaseTrait<AddForce<TForce...>>  {

    using Stencil = typename TTrait::Stencil;

    using Boundaries = typename TTrait::Boundaries;

    using PreProcessors = typename TTrait::PreProcessors;

    using PostProcessors = typename TTrait::PostProcessors;

    using Forces = tuple_cat_t<typename TTrait::Forces, std::tuple<TForce...>>;

    template<class TStencil>
    using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

    using Lattice = typename TTrait::Lattice;

    static constexpr int NumberOfComponents = TTrait::NumberOfComponents;

  };

  template<class... TForce>
  struct SetForce : BaseTrait<SetForce<TForce...>> {

    using Stencil = typename TTrait::Stencil;

    using Boundaries = typename TTrait::Boundaries;

    using PreProcessors = typename TTrait::PreProcessors;

    using PostProcessors = typename TTrait::PostProcessors;

    using Forces = std::tuple<TForce...>;

    template<class TStencil>
    using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

    using Lattice = typename TTrait::Lattice;

    static constexpr int NumberOfComponents = TTrait::NumberOfComponents;

  };

  template<class... TBoundary>
  struct AddBoundary : BaseTrait<AddBoundary<TBoundary...>> {

    using Stencil = typename TTrait::Stencil;

    using Boundaries = tuple_cat_t<typename TTrait::Boundaries, std::tuple<TBoundary...>>;

    using PreProcessors = typename TTrait::PreProcessors;

    using PostProcessors = typename TTrait::PostProcessors;

    using Forces = typename TTrait::Forces;
    
    template<class TStencil>
    using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

    using Lattice = typename TTrait::Lattice;

    static constexpr int NumberOfComponents = TTrait::NumberOfComponents;

  };

  template<class... TBoundary>
  struct SetBoundary : BaseTrait<SetBoundary<TBoundary...>> {

    using Stencil = typename TTrait::Stencil;

    using Boundaries = std::tuple<TBoundary...>;

    using PreProcessors = typename TTrait::PreProcessors;

    using PostProcessors = typename TTrait::PostProcessors;

    using Forces = typename TTrait::Forces;

    template<class TStencil>
    using CollisionModel = typename TTrait::template CollisionModel<TStencil>;

    using Lattice = typename TTrait::Lattice;

    static constexpr int NumberOfComponents = TTrait::NumberOfComponents;

  };

  template<class TStencil>
  struct SetStencil : BaseTrait<SetStencil<TStencil>> {

    using Stencil = TStencil;

    using Boundaries = typename TTrait::Boundaries;

    using PreProcessors = typename TTrait::PreProcessors;

    using PostProcessors = typename TTrait::PostProcessors;

    using Forces = typename TTrait::Forces;

    template<class stencil1>
    using CollisionModel = typename TTrait::template CollisionModel<stencil1>;

    using Lattice = typename TTrait::Lattice;

    static constexpr int NumberOfComponents = TTrait::NumberOfComponents;

  };

  template<template<class> class TModel>
  struct SetCollisionOperator : BaseTrait<SetCollisionOperator<TModel>> {

    using Stencil = typename TTrait::Stencil;

    using Boundaries = typename TTrait::Boundaries;

    using PreProcessors = typename TTrait::PreProcessors;

    using PostProcessors = typename TTrait::PostProcessors;

    using Forces = typename TTrait::Forces;

    template<class TStencil>
    using CollisionModel = TModel<TStencil>;

    using Lattice = typename TTrait::Lattice;

    static constexpr int NumberOfComponents = TTrait::NumberOfComponents;

  };

};

template<std::size_t i, class TTuple, std::size_t... is>
constexpr auto element_as_tuple(const TTuple tuple, std::index_sequence<is...>)
{
    if constexpr (!(std::is_same_v<std::tuple_element_t<i, TTuple>, 
                  std::tuple_element_t<is, TTuple>> || ...))
        return std::make_tuple(std::get<i>(tuple));
    else
        return std::make_tuple();
}

template<class TTuple, std::size_t... is>
constexpr auto make_tuple_unique(const TTuple tuple, std::index_sequence<is...>)
{
    return std::tuple_cat(element_as_tuple<is>(tuple, 
                          std::make_index_sequence<is>{})...);
}

template<class... Tuples>
constexpr auto make_tuple_unique(const Tuples... tuples)
{
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

struct Cartesian{};
struct AllDirections{};
struct One{};

template <typename TKey, typename TValue>
struct kv
{
    using Key = TKey;

    static TValue Value;
};

template <typename TKey, typename TValue>
TValue kv<TKey,TValue>::Value;

template <typename...>
struct ct_map;

template<>
struct ct_map<>
{

    template<typename>
    struct keyexists
    {
        static constexpr bool exists = false;
    };

    template<typename T>
    struct get
    {
        static inline constexpr int noKey(){
          if constexpr (!sizeof(T)){
            static_assert(!!sizeof(T), "Key does not exist in map.");
          }
          return 0;
        }
        
        static constexpr auto val=noKey();
    };
};

template<typename TKey, typename TValue, typename... TRest>
struct ct_map<kv<TKey, TValue>, TRest...>
{
    template<typename TKKey>
    struct keyexists
    {
        static constexpr bool exists = 
            (std::is_same<typename std::remove_reference<TKKey>::type, typename std::remove_reference<TKey>::type>::value) ?
            true :
            ct_map<TRest...>::template keyexists<TKKey>::exists;
    };

    template<typename TKKey>
    struct get
    {

        static inline constexpr auto& findVal(){
          static_assert(sizeof...(TRest)!=0||std::is_same<typename std::remove_reference<TKKey>::type, typename std::remove_reference<TKey>::type>::value, "Key does not exist in map.");
          if constexpr (sizeof...(TRest)!=0){
            return (std::is_same<TKKey, TKey>::value) ?
              kv<TKey, TValue>::Value : ct_map<TRest...>::template get<TKKey>::val;
          }
          else {
            return kv<TKey, TValue>::Value;
          }
          
        }

        static constexpr auto& val = findVal();
    };
};



template <typename TKey, typename TValue>
struct kv_types
{
  using Key = TKey;

  using Type = TValue;
};

template <typename...>
struct ct_map_types;

template<>
struct ct_map_types<>
{

    template<typename>
    struct keyexists
    {
        static constexpr bool exists = false;
    };

    template<typename T>
    struct get
    {
        static inline constexpr int noKey(){
          if constexpr (!sizeof(T)){
            static_assert(!!sizeof(T), "Key does not exist in map.");
          }
          return 0;
        }
        
        using TValue = decltype(noKey());

        TValue val = noKey();
    };
};

template<typename TKey, typename TValue, typename... TRest>
struct ct_map_types<kv<TKey, TValue>, TRest...>
{
    template<typename TKKey>
    struct keyexists
    {
        static constexpr bool exists = 
            (std::is_same<typename std::remove_reference<TKKey>::type, typename std::remove_reference<TKey>::type>::value) ?
            true :
            ct_map_types<TRest...>::template keyexists<TKKey>::exists;
    };

    template<typename TKKey>
    struct get
    {

        static inline constexpr auto findVal(){
          static_assert(sizeof...(TRest)!=0||std::is_same<typename std::remove_reference<TKKey>::type, typename std::remove_reference<TKey>::type>::value, "Key does not exist in map.");
          if constexpr (sizeof...(TRest)!=0){
            typename std::conditional_t<(std::is_same<TKKey, TKey>::value),typename kv_types<TKey, TValue>::Type,typename ct_map_types<TRest...>::template get<TKKey>::Type> val = {};
            return val;
          }
          else {
            typename kv_types<TKey, TValue>::Type val = {};
            return val;
          }
          
        }

        using Type = decltype(findVal());

        Type val = findVal();
    };
};


/*
template<typename ... input_t>
using tuple_cat_t=
decltype(std::tuple_cat(
    std::declval<input_t>()...
));

template<class TTrait, class TForce>
auto add_force(){return generate_trait<TTrait::Stencil,TTrait::Boundaries,tuple_cat_t::,TTrait::Properties>;}
*/


//BELOW: Not currently used implementations of type erasure, recursive template loops, MRT generation and a class
//instance counter.

/*
class LBModel {
	 
public:
    template <typename T>
    LBModel(T&& obj): object(std::make_shared<Model<T>>(obj)){}
      
    void collide() const {
        object->collide(); 
    }
    void precompute() const{
        object->precompute();
    }

    void initialise() const{
        object->initialise();
    }
	
    struct Concept {
        virtual ~Concept() {}
        virtual void collide() const = 0;
        virtual void precompute() const = 0;
        virtual void initialise() const = 0;
    };

   template<typename T>
   struct Model : Concept {
       Model(const T& t) : object(t) {}
	   void collide() const override {
		   object.collide();
	   }
       void precompute() const override {
		   object.precompute();
	   }
       void initialise() const override {
		   object.initialise();
	   }
     private:
       T object;
   };
private:
   std::shared_ptr<const Concept> object;
};

template<int index,typename T,typename... args>
struct TemplateLoop{
    template<void givenFunction(T* val,args... a)>
    inline static const void runGivenFunction(T* val,args... arguments){
        givenFunction<index>(val+1,(arguments,...));
        TemplateLoop<index-1,T,args...>::runGivenFunction<givenFunction>(val,(arguments,...));
    }
};

template<typename T,typename... args>
struct TemplateLoop<0,T,args...>{
    template<void givenFunction(T* val,args... arguments)>
    inline static const void runGivenFunction(T* val,args... arguments){
        
    }
};

template<typename T1,template<typename T2> class TStencil>
struct GenerateMRT{
    static const int mQ=TStencil<T1>::mQ;
    int mMRTMatrix[mQ][mQ];
    
    constexpr GenerateMRT():mMRTMatrix(){
        for (int i=0;i<mQ;i++){
            setMRTMatrixColumn<0>(i);
        }
    }

    template<int columnidx>
    inline constexpr void setMRTMatrixColumn(int rowidx){
        mMRTMatrix[rowidx][columnidx]=TStencil<T1>::template mMoments<columnidx>[rowidx];
        setMRTMatrixColumn<columnidx+1>(rowidx);
    }

    template<>
    inline constexpr void setMRTMatrixColumn<mQ>(int i){
        
    }
};

struct MRT{

};

template <typename T>
struct Counter
{
    Counter() {++counter;}
    virtual ~Counter() {--counter;}
    static int counter;
};
template <typename classtocount> int Counter<classtocount>::counter(0);
*/
