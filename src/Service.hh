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


template<class lattice>
inline int computeXGlobal(const int k) //Compute X direction from a given k, the convention in this code is that
                          //k will iterate over the z direction first, then increment y by 1 once it reaches LZ,
                          //then repeat the iteration over z. Once it reaches LY x will be incremented and this
                          //process continues
{  

  return lattice::LXMPIOffset+int((k-lattice::HaloSize)/(float) (lattice::LZ*lattice::LY));


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

template<class force>
typename force::Method getMethod(force& f){
    return std::declval<typename force::Method>();
}

template<class force>
decltype(force::Prefactor) getForcePrefactor(force& f){
    return std::declval<decltype(force::Prefactor)>();
}

void print() {
    if (mpi.rank != 0) return;
    std::cout << std::endl;
}

template <typename T, typename ... Args>
void print(std::vector<T> first, Args ... args) {
    if (mpi.rank != 0) return;
    for (auto elem : first) {
      std::cout << elem << " ";
    }
    print(args...);
}

template <typename T, typename ... Args>
void print(T first, Args ... args) {
    if (mpi.rank != 0) return;
    std::cout << first << " ";
    print(args...);
}

void printAll() {
    std::cout << std::endl;
}

template <typename T, typename ... Args>
void printAll(std::vector<T> first, Args ... args) {
    for (auto elem : first) {
      std::cout << elem << " ";
    }
    printAll(args...);
}

template <typename T, typename ... Args>
void printAll(T first, Args ... args) {
    std::cout << first << " ";
    printAll(args...);
}


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

template<typename wanted_type, typename T> struct is_wanted_type;

template<typename wanted_type, typename ...Types>
struct is_wanted_type<wanted_type, std::tuple<Types...>> {

    static constexpr bool wanted = (std::is_same_v<wanted_type, Types> || ...);

};

// Ok, the ith index in the tuple, here's its std::tuple_element type.
// And wanted_element_t is a tuple of all types we want to extract.
//
// Based on which way the wind blows we'll produce either a std::tuple<>
// or a std::tuple<tuple_element_t>.

template<size_t i, typename tuple_element_t,
         typename wanted_element_t,
         bool wanted = is_wanted_type<tuple_element_t, wanted_element_t>::wanted>
struct extract_type {

    template<typename tuple_type>
    inline static auto do_extract_type(tuple_type &t) {

        return std::tuple<>{};

    }

};


template<size_t i, typename tuple_element_t, typename wanted_element_t>
struct extract_type<i, tuple_element_t, wanted_element_t, true> {

    template<typename tuple_type>
    inline static auto do_extract_type(tuple_type &t) {

        return std::tie(std::get<i>(t));

    }

};

// And now, a simple fold expression to pull out all wanted types
// and tuple-cat them together.

template<typename wanted_element_t, typename tuple_type, size_t ...i>
inline auto get_type_t(tuple_type &t, std::index_sequence<i...>) {

    return std::tuple_cat(extract_type<i,
                          typename std::tuple_element<i, tuple_type>::type,
                          wanted_element_t>::do_extract_type(t)...);

}


template<typename ...wanted_element_t, typename ...types>
inline auto get_type(std::tuple<types...> &t) {

    return get_type_t<std::tuple<wanted_element_t...>>(t, std::make_index_sequence<sizeof...(types)>());
        
}

template<class Base, typename Tuple>
struct CheckBase;

template<class Base, typename... Types>
struct CheckBase<Base, std::tuple<Types...>> : std::conjunction<std::is_base_of<Base,Types>...> {};


template <template <typename...> class C, typename...Ts>
std::true_type is_base_of_template_impl(const C<Ts...>*);

template <template <typename...> class C>
std::false_type is_base_of_template_impl(...);

template <template <typename...> class C,typename T>
using is_base_of_template = decltype(is_base_of_template_impl<C>(std::declval<T*>()));

template<template<class> class Base, typename Tuple>
struct CheckBaseTemplate;

template<template<class> class Base, typename... Types>
struct CheckBaseTemplate<Base, std::tuple<Types...>> : std::conjunction<is_base_of_template<Base,Types>...> {};

template<typename ... input_t>
using tuple_cat_t=
decltype(std::tuple_cat(
    std::declval<input_t>()...
));

template<class trait>
struct BaseTrait{

  template<class... preprocessor>
  struct AddPreProcessor : BaseTrait<AddPreProcessor<preprocessor...>>  {

    using Stencil = typename trait::Stencil;

    using Boundaries = typename trait::Boundaries;

    using PreProcessors = tuple_cat_t<typename trait::PreProcessors, std::tuple<preprocessor...>>;

    using PostProcessors = typename trait::PostProcessors;

    using Forces = typename trait::Forces;

    template<class stencil>
    using CollisionModel = typename trait::template CollisionModel<stencil>;

    using Lattice = typename trait::Lattice;

    static constexpr int NumberOfComponents = trait::NumberOfComponents;

  };

  template<class... preprocessor>
  struct SetPreProcessor : BaseTrait<SetPreProcessor<preprocessor...>> {

    using Stencil = typename trait::Stencil;

    using Boundaries = typename trait::Boundaries;

    using PreProcessors = std::tuple<preprocessor...>;

    using PostProcessors = typename trait::PostProcessors;

    using Forces = typename trait::Forces;

    template<class stencil>
    using CollisionModel = typename trait::template CollisionModel<stencil>;

    using Lattice = typename trait::Lattice;

    static constexpr int NumberOfComponents = trait::NumberOfComponents;

  };

  template<class... postprocessor>
  struct AddPostProcessor : BaseTrait<AddPostProcessor<postprocessor...>>  {

    using Stencil = typename trait::Stencil;

    using Boundaries = typename trait::Boundaries;

    using PreProcessors = typename trait::PreProcessors;

    using PostProcessors = tuple_cat_t<typename trait::PostProcessors, std::tuple<postprocessor...>>;

    using Forces = typename trait::Forces;

    template<class stencil>
    using CollisionModel = typename trait::template CollisionModel<stencil>;

    using Lattice = typename trait::Lattice;

    static constexpr int NumberOfComponents = trait::NumberOfComponents;

  };

  template<class... postprocessor>
  struct SetPostProcessor : BaseTrait<SetPostProcessor<postprocessor...>> {

    using Stencil = typename trait::Stencil;

    using Boundaries = typename trait::Boundaries;

    using PreProcessors = typename trait::PreProcessors;

    using PostProcessors = std::tuple<postprocessor...>;

    using Forces = typename trait::Forces;

    template<class stencil>
    using CollisionModel = typename trait::template CollisionModel<stencil>;

    using Lattice = typename trait::Lattice;

    static constexpr int NumberOfComponents = trait::NumberOfComponents;

  };
  
  template<class... force>
  struct AddForce : BaseTrait<AddForce<force...>>  {

    using Stencil = typename trait::Stencil;

    using Boundaries = typename trait::Boundaries;

    using PreProcessors = typename trait::PreProcessors;

    using PostProcessors = typename trait::PostProcessors;

    using Forces = tuple_cat_t<typename trait::Forces, std::tuple<force...>>;

    template<class stencil>
    using CollisionModel = typename trait::template CollisionModel<stencil>;

    using Lattice = typename trait::Lattice;

    static constexpr int NumberOfComponents = trait::NumberOfComponents;

  };

  template<class... force>
  struct SetForce : BaseTrait<SetForce<force...>> {

    using Stencil = typename trait::Stencil;

    using Boundaries = typename trait::Boundaries;

    using PreProcessors = typename trait::PreProcessors;

    using PostProcessors = typename trait::PostProcessors;

    using Forces = std::tuple<force...>;

    template<class stencil>
    using CollisionModel = typename trait::template CollisionModel<stencil>;

    using Lattice = typename trait::Lattice;

    static constexpr int NumberOfComponents = trait::NumberOfComponents;

  };

  template<class... boundary>
  struct AddBoundary : BaseTrait<AddBoundary<boundary...>> {

    using Stencil = typename trait::Stencil;

    using Boundaries = tuple_cat_t<typename trait::Boundaries, std::tuple<boundary...>>;

    using PreProcessors = typename trait::PreProcessors;

    using PostProcessors = typename trait::PostProcessors;

    using Forces = typename trait::Forces;
    
    template<class stencil>
    using CollisionModel = typename trait::template CollisionModel<stencil>;

    using Lattice = typename trait::Lattice;

    static constexpr int NumberOfComponents = trait::NumberOfComponents;

  };

  template<class... boundary>
  struct SetBoundary : BaseTrait<SetBoundary<boundary...>> {

    using Stencil = typename trait::Stencil;

    using Boundaries = std::tuple<boundary...>;

    using PreProcessors = typename trait::PreProcessors;

    using PostProcessors = typename trait::PostProcessors;

    using Forces = typename trait::Forces;

    template<class stencil>
    using CollisionModel = typename trait::template CollisionModel<stencil>;

    using Lattice = typename trait::Lattice;

    static constexpr int NumberOfComponents = trait::NumberOfComponents;

  };

  template<class stencil>
  struct SetStencil : BaseTrait<SetStencil<stencil>> {

    using Stencil = stencil;

    using Boundaries = typename trait::Boundaries;

    using PreProcessors = typename trait::PreProcessors;

    using PostProcessors = typename trait::PostProcessors;

    using Forces = typename trait::Forces;

    template<class stencil1>
    using CollisionModel = typename trait::template CollisionModel<stencil1>;

    using Lattice = typename trait::Lattice;

    static constexpr int NumberOfComponents = trait::NumberOfComponents;

  };

  template<template<class> class model>
  struct SetCollisionModel : BaseTrait<SetCollisionModel<model>> {

    using Stencil = typename trait::Stencil;

    using Boundaries = typename trait::Boundaries;

    using PreProcessors = typename trait::PreProcessors;

    using PostProcessors = typename trait::PostProcessors;

    using Forces = typename trait::Forces;

    template<class stencil>
    using CollisionModel = model<stencil>;

    using Lattice = typename trait::Lattice;

    static constexpr int NumberOfComponents = trait::NumberOfComponents;

  };

};

template<std::size_t i, class Tuple, std::size_t... is>
constexpr auto element_as_tuple(const Tuple tuple, std::index_sequence<is...>)
{
    if constexpr (!(std::is_same_v<std::tuple_element_t<i, Tuple>, 
                  std::tuple_element_t<is, Tuple>> || ...))
        return std::make_tuple(std::get<i>(tuple));
    else
        return std::make_tuple();
}

template<class Tuple, std::size_t... is>
constexpr auto make_tuple_unique(const Tuple tuple, std::index_sequence<is...>)
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



template <typename T, typename Tuple>
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

template <typename key, typename valuetype>
struct kv
{
    using Key = key;

    static valuetype Value;
};

template <typename key, typename valuetype>
valuetype kv<key,valuetype>::Value;

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

template<typename key, typename value, typename... rest>
struct ct_map<kv<key, value>, rest...>
{
    template<typename kkey>
    struct keyexists
    {
        static constexpr bool exists = 
            (std::is_same<typename std::remove_reference<kkey>::type, typename std::remove_reference<key>::type>::value) ?
            true :
            ct_map<rest...>::template keyexists<kkey>::exists;
    };

    template<typename kkey>
    struct get
    {

        static inline constexpr auto& findVal(){
          static_assert(sizeof...(rest)!=0||std::is_same<typename std::remove_reference<kkey>::type, typename std::remove_reference<key>::type>::value, "Key does not exist in map.");
          if constexpr (sizeof...(rest)!=0){
            return (std::is_same<kkey, key>::value) ?
              kv<key, value>::Value : ct_map<rest...>::template get<kkey>::val;
          }
          else {
            return kv<key, value>::Value;
          }
          
        }

        static auto& val = findVal();
    };
};



template <typename key, typename valuetype>
struct kv_types
{
  using Key = key;

  using Type = valuetype;
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
        
        using valuetype = decltype(noKey());

        valuetype val = noKey();
    };
};

template<typename key, typename value, typename... rest>
struct ct_map_types<kv<key, value>, rest...>
{
    template<typename kkey>
    struct keyexists
    {
        static constexpr bool exists = 
            (std::is_same<typename std::remove_reference<kkey>::type, typename std::remove_reference<key>::type>::value) ?
            true :
            ct_map_types<rest...>::template keyexists<kkey>::exists;
    };

    template<typename kkey>
    struct get
    {

        static inline constexpr auto findVal(){
          static_assert(sizeof...(rest)!=0||std::is_same<typename std::remove_reference<kkey>::type, typename std::remove_reference<key>::type>::value, "Key does not exist in map.");
          if constexpr (sizeof...(rest)!=0){
            typename std::conditional_t<(std::is_same<kkey, key>::value),typename kv_types<key, value>::Type,typename ct_map_types<rest...>::template get<kkey>::Type> val = {};
            return val;
          }
          else {
            typename kv_types<key, value>::Type val = {};
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

template<class trait, class force>
auto add_force(){return generate_trait<trait::Stencil,trait::Boundaries,tuple_cat_t::,trait::Properties>;}
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

template<typename T1,template<typename T2> class stencil>
struct GenerateMRT{
    static const int m_Q=stencil<T1>::m_Q;
    int m_MRTMatrix[m_Q][m_Q];
    
    constexpr GenerateMRT():m_MRTMatrix(){
        for (int i=0;i<m_Q;i++){
            setMRTMatrixColumn<0>(i);
        }
    }

    template<int columnidx>
    inline constexpr void setMRTMatrixColumn(int rowidx){
        m_MRTMatrix[rowidx][columnidx]=stencil<T1>::template m_Moments<columnidx>[rowidx];
        setMRTMatrixColumn<columnidx+1>(rowidx);
    }

    template<>
    inline constexpr void setMRTMatrixColumn<m_Q>(int i){
        
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
