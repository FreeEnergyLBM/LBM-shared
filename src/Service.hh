#ifndef SERVICE_HEADER
#define SERVICE_HEADER

#include <map>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <any>
#include<charconv>
#include <memory>

template<class model>
void LBMPrint(model& m){//THIS IS TEMPORARY
    std::cout<<"Distributions: "<<std::flush;
    for(int i=0;i<9;i++)std::cout<<m.getDistribution()[i+9*15]<<" "<<std::flush;
    std::cout<<std::endl;
    std::cout<<"Density: "<<m.getDensity(15)<<" "<<std::endl;
    std::cout<<"Velocity: "<<std::flush;
    for(int i=0;i<2;i++)std::cout<<m.getVelocity()[15*2+i]<<" "<<std::flush;
    std::cout<<std::endl;
}

int computeX(const int k)
{  

  return int(k/(float) (LZ*LY));


}

int computeY(const int k)
{  

  return int((k-computeX(k)*LZ*LY)/(float) LZ);

}

int computeZ(const int k)
{  

  return k-computeX(k)*LZ*LY-computeY(k)*LZ;

}


template<class stencil,template <typename givenstencil,int id> class data,class... forces>
struct trait{
    using Stencil=stencil;  
    template<int id>
    using Data=data<Stencil,id> ;
    typedef std::tuple<forces&...> Forces;
};


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