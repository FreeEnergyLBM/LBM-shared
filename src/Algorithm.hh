#ifndef ALGORITHM_HEADER
#define ALGORITHM_HEADER
#include <tuple>

template<class ...Model>
class Algorithm{
    public:
        Algorithm(Model&... Models):mt_Models(Models...){

        }
        void evolve();
        void initialise();
    private:
        void precomputeStep();
        void calculateCollisionStep();
        void calculateMomentaStep();
        std::tuple<Model&...> mt_Models;
};

template<class ...Model>
void Algorithm<Model...>::precomputeStep(){
    if constexpr (sizeof...(Model)!=0) (std::get<Model&>(mt_Models).precompute(),...);
    else;
}

template<class ...Model>
void Algorithm<Model...>::calculateCollisionStep(){
    if constexpr (sizeof...(Model)!=0) (std::get<Model&>(mt_Models).collide(),...);
    else;
}

template<class ...Model>
void Algorithm<Model...>::calculateMomentaStep(){
    if constexpr (sizeof...(Model)!=0) (std::get<Model&>(mt_Models).computeMomenta(),...);
    else;
}

template<class ...Model>
void Algorithm<Model...>::evolve(){
    precomputeStep();
    calculateCollisionStep();
    calculateMomentaStep();
}

template<class ...Model>
void Algorithm<Model...>::initialise(){
    if constexpr (sizeof...(Model)!=0) (std::get<Model&>(mt_Models).initialise(),...);
    else;
}

#endif