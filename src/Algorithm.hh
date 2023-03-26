#ifndef ALGORITHM_HEADER
#define ALGORITHM_HEADER
#include <tuple>
#include <iostream>

//Algorithm.hh: This file contains the Algorithm class, which will take the given models and provide functions
//to initialise the models and evolve the LBM algorithm with the models, stencils, data types, forcing terms
//and boundary conditions selected.

template<class ...Model> //Pass any number of models as to this class (the "..." allows this)
class Algorithm{
    public:

        Algorithm(Model&... Models):mt_Models(Models...){} //Pass objects of each model to this class when it is
                                                           //constructed. The order of objects passed must
                                                           //correspond to the order of models pased in the
                                                           //initial template.

        void evolve(); //Evolve the LBM algorithm. For now this is one timestep

        void initialise(); //Perform necessary initialisation

    private:

        void precomputeStep(); //Perform any necessary calculations before collision can take place at each
                               //timestep
        
        void calculateCollisionStep(); //Calculate collision for each model

        void calculateBoundaryStep(); //Apply boundary conditions for each model

        void calculateMomentaStep(); //Calculate momenta (density, velocity) for each model.

        std::tuple<Model&...> mt_Models; //Tuple containing objects of each model

};

template<class ...Model>
void Algorithm<Model...>::evolve(){ //See above

    precomputeStep();

    calculateCollisionStep();

    calculateBoundaryStep();

    calculateMomentaStep();
    
}

template<class ...Model>
void Algorithm<Model...>::precomputeStep(){

    if constexpr (sizeof...(Model)!=0){ //Check if we have any models in the class
        std::apply([](Model&... models){ //std::apply will run whatever is contained within the {} by expanding
                                         //the tuple and using the elements as arguments. This is neccessary
                                         //as std::get<> does not work if we have repeats in "Model&...".
                (models.precompute(),...); //This expression runs models.precompute() for every "Model&..."
            }, mt_Models);
    }
    else; //THROW SOME ERROR, ALSO PROBABLY JUST DO THIS ONCE

}

template<class ...Model>
void Algorithm<Model...>::calculateCollisionStep(){ //...
    
    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... models){
                (models.collide(),...);
            }, mt_Models);
    }
    else;

}

template<class ...Model>
void Algorithm<Model...>::calculateBoundaryStep(){ //...
    
    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... models){
                (models.boundaries(),...);
            }, mt_Models);
    }
    else;

}

template<class ...Model>
void Algorithm<Model...>::calculateMomentaStep(){ //...

    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... models){
                (models.computeMomenta(),...);
            }, mt_Models);
    }
    else;

}

template<class ...Model>
void Algorithm<Model...>::initialise(){ //...

    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... models){
                (models.initialise(),...);
            }, mt_Models);
    }
    else;

}

#endif