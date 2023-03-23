#ifndef ALGORITHM_HEADER
#define ALGORITHM_HEADER
#include <tuple>
#include <iostream>

template<class ...Model>
class Algorithm{
    public:

        Algorithm(Model&... Models):mt_Models(Models...){}

        void evolve();

        void initialise();

    private:

        void precomputeStep();
        
        void calculateCollisionStep();

        void calculateBoundaryStep();

        void calculateMomentaStep();

        std::tuple<Model&...> mt_Models; //CHANGE THIS TO TUPLE OF VECTORS

};

template<class ...Model>
void Algorithm<Model...>::evolve(){

    precomputeStep();

    calculateCollisionStep();

    calculateMomentaStep();
    
}

template<class ...Model>
void Algorithm<Model...>::precomputeStep(){

    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... tests){
                (tests.precompute(),...);
            }, mt_Models);
    }
    else;

}

template<class ...Model>
void Algorithm<Model...>::calculateCollisionStep(){
    
    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... tests){
                (tests.collide(),...);
            }, mt_Models);
    }
    else;

}

template<class ...Model>
void Algorithm<Model...>::calculateBoundaryStep(){
    
    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... tests){
                (tests.boundaries(),...);
            }, mt_Models);
    }
    else;

}

template<class ...Model>
void Algorithm<Model...>::calculateMomentaStep(){

    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... tests){
                (tests.computeMomenta(),...);
            }, mt_Models);
    }
    else;

}

template<class ...Model>
void Algorithm<Model...>::initialise(){

    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... tests){
                (tests.initialise(),...);
            }, mt_Models);
    }
    else;

}

#endif