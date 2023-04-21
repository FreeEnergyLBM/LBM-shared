#ifndef ALGORITHM_HEADER
#define ALGORITHM_HEADER
#include <tuple>
#include <iostream>

/**
 * \file Algorithm.hh
 * \brief This file contains classes to run the lattice boltzmann algorithm for the chosen models.
 * The classes in this file will take the given models and provide functions
 * to initialise the models and evolve the LBM algorithm with the models, stencils, data types, forcing terms
 * and boundary conditions selected. If a different algorithm is required, add a new class here and inherit from
 * the class Algorithm.
 */

/**
 * \brief This class takes runs the standard LBM algorithm for each of the models specified in the template.
 * The Algorithm class takes any number of model classes as template arguments. It will put these in a tuple
 * (mt_Models) and then use the std::apply function to perform the compute, boundary, momenta and precompute
 * calculations for each model. The models that are passed to the template must therefore have these public
 * functions. The class must be given objects of each model when it is constructed. Then, the models can be
 * initialised and evolved by one timestep at a time.
 */
template<class ...Model>
class Algorithm{
    public:

        /**
         * \brief Constructor for the class that will fill the tuple "mt_Models" with given objects of each model.
         * This constructor will accept model objects corresponding to each model in the order they are given in the
         * template arguments. This is then used to initialise the tuple "mt_Models" with references to each
         * object. Note that models will be computed in the order they are specified.
         * \param Models Objects of each model in the order specified by the template parameter "...Model".
         */
        Algorithm(Model&... Models):mt_Models(Models...){}

        /**
         * \brief Constructor for the class that will fill the tuple "mt_Models" with given objects of each model.
         * This constructor will accept model objects corresponding to each model in the order they are given in the
         * template arguments. This is then used to initialise the tuple "mt_Models" with references to each
         * object. Note that models will be computed in the order they are specified.
         * \param Models Objects of each model in the order specified by the template parameter "...Model".
         */
        Algorithm():mt_Models(*new Model...){}

        /**
         * \brief Function that will evolve the lattice Boltzmann algorithm by one timestep.
         */
        void evolve();

        /**
         * \brief Function that will perform necessary initialisations for each model (e.g. set distributions to equilibrium).
         */
        void initialise();

    private:

        /**
         * \brief Perform any necessary calculations before collision can take place at each timestep.
         */
        void precomputeStep(); 
        
        /**
         * \brief Calculate collision (and streaming currently) for each model over the entire lattice.
         */
        void calculateCollisionStep(); 

        /**
         * \brief Apply boundary conditions for each model over the entire lattice.
         */
        void calculateBoundaryStep();

        /**
         * \brief Calculate momenta (density, velocity) for each model over the entire lattice.
         */
        void calculateMomentaStep();

        /*
         * Tuple containing references to objects of each Model... passed through the constructor.
         */
        std::tuple<Model&...> mt_Models;

};

/**
 * \details This function will perform the precompute step (e.g. gradient calculations, chemical potential
 *          calculations) then the collision (and streaming currently) step, then the boundary calculations
 *          and finally it will compute the macroscopic variables (density, velocity etc.). It will do this
 *          for every model.
 */
template<class ...Model>
void Algorithm<Model...>::evolve(){

    precomputeStep();
    
    calculateCollisionStep();
    
    calculateBoundaryStep();
    
    calculateMomentaStep();
    
}

/**
 * \details This function first checks if the algoritm has any models in the template arguments. Then, the 
 *          std::apply() function is used to apply a lambda function, taking arguments as the models stored in the
 *          tuple "mt_Models". In this case, the lambda function applies "(models.initialise(),...);", which will 
 *          run the "initialise()" function for every model in the tuple. This function might set distributions to
 *          equilibrium and set macroscopic variables to some initial value, for instance.
 */
template<class ...Model>
void Algorithm<Model...>::initialise(){ //...

    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... models){
                (models.initialise(),...);
            }, mt_Models);
    }
    else;

}

/**
 * \details This function first checks if the algoritm has any models in the template arguments. Then, the 
 *          std::apply() function is used to apply a lambda function, taking arguments as the models stored in the
 *          tuple "mt_Models". In this case, the lambda function applies "(models.precompute(),...);", which will 
 *          run the "precompute()" function for every model in the tuple. This function might perform some gradient
 *          calculations needed in the forcing terms, for instance.
 */
template<class ...Model>
void Algorithm<Model...>::precomputeStep(){

    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... models){
                (models.precompute(),...);
            }, mt_Models);
    }
    else;

}

/**
 * \details This function first checks if the algoritm has any models in the template arguments. Then, the 
 *          std::apply() function is used to apply a lambda function, taking arguments as the models stored in the
 *          tuple "mt_Models". In this case, the lambda function applies "(models.collide(),...);", which will 
 *          run the "collide()" function for every model in the tuple. This function will collide based on the 
 *          chosen collision model and will also perform streaming.
 */
template<class ...Model>
void Algorithm<Model...>::calculateCollisionStep(){ //...
    
    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... models){
                (models.collide(),...);
            }, mt_Models);
    }
    else;

}

/**
 * \details This function first checks if the algoritm has any models in the template arguments. Then, the 
 *          std::apply() function is used to apply a lambda function, taking arguments as the models stored in the
 *          tuple "mt_Models". In this case, the lambda function applies "(models.boundaries(),...);", which will 
 *          run the "boundaries()" function for every model in the tuple. This function might apply bounceback and
 *          outflow boundaries, depending on the geometry labels, for instance.
 */
template<class ...Model>
void Algorithm<Model...>::calculateBoundaryStep(){ //...
    
    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... models){
                (models.boundaries(),...);
            }, mt_Models);
    }
    else;

}

/**
 * \details This function first checks if the algoritm has any models in the template arguments. Then, the 
 *          std::apply() function is used to apply a lambda function, taking arguments as the models stored in the
 *          tuple "mt_Models". In this case, the lambda function applies "(models.computeMomenta(),...);", which
 *          will run the "computeMomenta()" function for every model in the tuple. This function might set
 *          calculate density and velocity based on the distributions, for instance.
 */
template<class ...Model>
void Algorithm<Model...>::calculateMomentaStep(){ //...

    if constexpr (sizeof...(Model)!=0){
        std::apply([](Model&... models){
                (models.computeMomenta(),...);
            }, mt_Models);
    }
    else;

}

#endif