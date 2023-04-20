#ifndef SERVICE_HEADER
#define SERVICE_HEADER

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
#ifdef MPIPARALLEL
#include <mpi.h>
#endif
//Sertive.hh: This will contain some commonly used functions with various uses.

int computeX(const int k) //Compute X direction from a given k, the convention in this code is that
                          //k will iterate over the z direction first, then increment y by 1 once it reaches LZ,
                          //then repeat the iteration over z. Once it reaches LY x will be incremented and this
                          //process continues
{  

  return int(k/(float) (LZ*LY));


}

int computeY(const int k) //Compute Y direction from a given k, this uses information from the X direction
{  

  return int((k-computeX(k)*LZ*LY)/(float) LZ);

}

int computeZ(const int k) //Compute Y direction from a given k, this uses information from the X and Y directions
{  

  return k-computeX(k)*LZ*LY-computeY(k)*LZ;

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
  } else if constexpr (std::is_same_v<T, signed char>) {
    mpi_type = MPI_SIGNED_CHAR;
  } else if constexpr (std::is_same_v<T, unsigned char>) {
    mpi_type = MPI_UNSIGNED_CHAR;
  } else if constexpr (std::is_same_v<T, wchar_t>) {
    mpi_type = MPI_WCHAR;
  } else if constexpr (std::is_same_v<T, signed short>) {
    mpi_type = MPI_SHORT;
  } else if constexpr (std::is_same_v<T, unsigned short>) {
    mpi_type = MPI_UNSIGNED_SHORT;
  } else if constexpr (std::is_same_v<T, signed int>) {
    mpi_type = MPI_INT;
  } else if constexpr (std::is_same_v<T, unsigned int>) {
    mpi_type = MPI_UNSIGNED;
  } else if constexpr (std::is_same_v<T, signed long int>) {
     mpi_type = MPI_LONG;
  } else if constexpr (std::is_same_v<T, unsigned long int>) {
    mpi_type = MPI_UNSIGNED_LONG;
  } else if constexpr (std::is_same_v<T, signed long long int>) {
    mpi_type = MPI_LONG_LONG;
  } else if constexpr (std::is_same_v<T, unsigned long long int>) {
    mpi_type = MPI_UNSIGNED_LONG_LONG;
  } else if constexpr (std::is_same_v<T, float>) {
    mpi_type = MPI_FLOAT;
  } else if constexpr (std::is_same_v<T, double>) {
    mpi_type = MPI_DOUBLE;
  } else if constexpr (std::is_same_v<T, long double>) {
    mpi_type = MPI_LONG_DOUBLE;
  } else if constexpr (std::is_same_v<T, std::int8_t>) {
    mpi_type = MPI_INT8_T;
  } else if constexpr (std::is_same_v<T, std::int16_t>) {
    mpi_type = MPI_INT16_T;
  } else if constexpr (std::is_same_v<T, std::int32_t>) {
    mpi_type = MPI_INT32_T;
  } else if constexpr (std::is_same_v<T, std::int64_t>) {
    mpi_type = MPI_INT64_T;
  } else if constexpr (std::is_same_v<T, std::uint8_t>) {
    mpi_type = MPI_UINT8_T;
  } else if constexpr (std::is_same_v<T, std::uint16_t>) {
    mpi_type = MPI_UINT16_T;
  } else if constexpr (std::is_same_v<T, std::uint32_t>) {
    mpi_type = MPI_UINT32_T;
  } else if constexpr (std::is_same_v<T, std::uint64_t>) {
    mpi_type = MPI_UINT64_T;
  } else if constexpr (std::is_same_v<T, bool>) {
    mpi_type = MPI_C_BOOL;
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    mpi_type = MPI_C_COMPLEX;
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    mpi_type = MPI_C_DOUBLE_COMPLEX;
  } else if constexpr (std::is_same_v<T, std::complex<long double>>) {
    mpi_type = MPI_C_LONG_DOUBLE_COMPLEX;
  }
	
  assert(mpi_type != MPI_DATATYPE_NULL);

  return mpi_type;    
}
#endif
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

   template< typename T >
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
#endif