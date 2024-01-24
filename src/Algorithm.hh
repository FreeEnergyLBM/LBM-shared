#pragma once
#include <tuple>
#include <iostream>
#include "Lattice.hh"
#include "LBModels/ModelBase.hh"

/**
 * \file Algorithm.hh
 * \brief This file contains classes to run the lattice boltzmann algorithm for
 * the chosen models. The classes in this file will take the given models and 
 * provide functions to initialise the models and evolve the LBM algorithm with
 * the models, stencils, data types, forcing terms and boundary conditions
 * selected. If a different algorithm is required, add a new class here and
 * inherit from the class Algorithm.
 */

/**
 * \brief This class takes runs the standard LBM algorithm for each of the 
 * models specified in the template. The Algorithm class takes any number of
 * model classes as template arguments. It will put these in a tuple
 * (mt_Models) and then use the std::apply function to perform the compute,
 * boundary, momenta and precompute calculations for each model. The models
 * that are passed to the template must therefore have these public functions.
 * The class must be given objects of each model when it is constructed. Then,
 * the models can be initialised and evolved by one timestep at a time.
 * \tparam TModel Any number of LBM model classes.
 */
template<class ...TModel>
class Algorithm {
    static_assert(std::conjunction<is_base_of_template<ModelBase,TModel>...>::value,
                  "ERROR: At least one LBM model chosen is not a model class.");
    public:

        /**
         * \brief Constructor for the class that will fill the tuple "mt_Models" with given objects of each model.
         * This constructor will accept model objects corresponding to each model in the order they are given in the
         * template arguments. This is then used to initialise the tuple "mt_Models" with references to each
         * object. Note that models will be computed in the order they are specified. This constructor also runs the
         * initialise() function, which will perform necessary initialisation for each model.
         * \param Models Objects of each model in the order specified by the template parameter "...TModel".
         */
        Algorithm(TModel&... Models) 
            : mt_Models(Models...)
        {initialise();}

        /**
         * \brief Constructor for the class that will fill the tuple "mt_Models" with objects of each model.
         * This constructor will create new model objects in the order they are given in the
         * template arguments. This is then used to initialise the tuple "mt_Models". Note that models will be computed 
         * in the order they are specified. This constructor also runs the
         * initialise() function, which will perform necessary initialisation for each model.
         * \param Models Objects of each model in the order specified by the template parameter "...TModel".
         */
        
        Algorithm()
            : mt_Models(*new TModel...)
        {initialise();}

        /**
         * \brief Function that will evolve the lattice Boltzmann algorithm by one timestep.
         */
        inline void evolve();

        /**
         * \brief Function that will perform necessary initialisations for each model (e.g. set distributions to equilibrium).
         */
        inline void initialise();

    private:

        /**
         * \brief Perform any necessary calculations before collision can take place at each timestep.
         */
        inline void processorStep(); 
        
        /**
         * \brief Calculate collision (and streaming currently) for each model over the entire lattice.
         */
        inline void calculateCollisionStep(); 

        inline void calculateStreamingStep(); 

        /**
         * \brief Apply boundary conditions for each model over the entire lattice.
         */
        inline void calculateBoundaryStep();

        /**
         * \brief Calculate momenta (density, velocity) for each model over the entire lattice.
         */
        inline void calculateMomentaStep();

        /**
         * \brief Postprocess.
         *//*
        inline void postprocessStep();*/

        /*
         * Tuple containing references to objects of each TModel... passed through the constructor.
         */
        std::tuple<TModel&...> mt_Models;

};

/**
 * \details This function will perform the precompute step (e.g. gradient calculations, chemical potential
 *          calculations) then the collision (and streaming currently) step, then the boundary calculations,
 *          then it will compute the macroscopic variables (density, velocity etc.) and finally it will do 
 *          the postprocessing step (which can do the same things as the precompute step). It will do this 
 *          for every model.
 */
template<class ...TModel>
inline void Algorithm<TModel...>::evolve() {
    
    #pragma omp parallel
    {   
    
    calculateCollisionStep();

    calculateStreamingStep();
    
    calculateBoundaryStep();
    
    calculateMomentaStep();
    
    processorStep();
    
    }
}

/**
 * \details This function first checks if the algoritm has any models in the template arguments. Then, the 
 *          std::apply() function is used to apply a lambda function, taking arguments as the models stored in the
 *          tuple "mt_Models". This is necessary to access the elements of a tuple where there may be duplicate types.
 *          In this case, the lambda function applies "(models.initialise(),...);", which will 
 *          run the "initialise()" function for every model in the tuple. This function might set distributions to
 *          equilibrium and set macroscopic variables to some initial value, for instance.
 */
template<class ...TModel>
inline void Algorithm<TModel...>::initialise() { //...

    if constexpr (sizeof...(TModel) != 0) {

        std::apply([](TModel&... models) {

            (models.initialise(), ...);

        }, mt_Models);

    }
    
    processorStep();

}

/**
 * \details This function first checks if the algoritm has any models in the template arguments. Then, the 
 *          std::apply() function is used to apply a lambda function, taking arguments as the models stored in the
 *          tuple "mt_Models". In this case, the lambda function applies "(models.precompute(),...);", which will 
 *          run the "precompute()" function for every model in the tuple. This function might perform some gradient
 *          calculations needed in the forcing terms, for instance.
 */
template<class ...TModel>
inline void Algorithm<TModel...>::processorStep() {

    if constexpr (sizeof...(TModel) != 0) {

        std::apply([](TModel&... models) {
            
            (models.computeProcessors(), ...);

        }, mt_Models);

    }

}

/**
 * \details This function first checks if the algoritm has any models in the template arguments. Then, the 
 *          std::apply() function is used to apply a lambda function, taking arguments as the models stored in the
 *          tuple "mt_Models". In this case, the lambda function applies "(models.postprocess(),...);", which will 
 *          run the "postprocess()" function for every model in the tuple. This function might perform some gradient
 *          calculations needed in the forcing terms, for instance.
 *//*
template<class ...TModel>
inline void Algorithm<TModel...>::postprocessStep() {

    if constexpr (sizeof...(TModel) != 0) {

        std::apply([](TModel&... models) {

            (models.postprocess(), ...);

        }, mt_Models);

    }

}
*/
/**
 * \details This function first checks if the algoritm has any models in the template arguments. Then, the 
 *          std::apply() function is used to apply a lambda function, taking arguments as the models stored in the
 *          tuple "mt_Models". In this case, the lambda function applies "(models.collide(),...);", which will 
 *          run the "collide()" function for every model in the tuple. This function will collide based on the 
 *          chosen collision model and will also perform streaming.
 */
template<class ...TModel>
inline void Algorithm<TModel...>::calculateCollisionStep() { //...
    
    if constexpr (sizeof...(TModel) != 0){

        std::apply([](TModel&... models) {

            (models.collide(), ...);

        }, mt_Models);

    }

}

template<class ...TModel>
inline void Algorithm<TModel...>::calculateStreamingStep() { //...
    
    if constexpr (sizeof...(TModel) != 0){

        std::apply([](TModel&... models) {

            (models.stream(), ...);

        }, mt_Models);

    }

}

/**
 * \details This function first checks if the algoritm has any models in the template arguments. Then, the 
 *          std::apply() function is used to apply a lambda function, taking arguments as the models stored in the
 *          tuple "mt_Models". In this case, the lambda function applies "(models.boundaries(),...);", which will 
 *          run the "boundaries()" function for every model in the tuple. This function might apply bounceback and
 *          outflow boundaries, depending on the geometry labels, for instance.
 */
template<class ...TModel>
inline void Algorithm<TModel...>::calculateBoundaryStep() { //...
    
    if constexpr (sizeof...(TModel) != 0){

        std::apply([](TModel&... models) {

            (models.boundaries(), ...);

        }, mt_Models);

    }

}

/**
 * \details This function first checks if the algoritm has any models in the template arguments. Then, the 
 *          std::apply() function is used to apply a lambda function, taking arguments as the models stored in the
 *          tuple "mt_Models". In this case, the lambda function applies "(models.computeMomenta(),...);", which
 *          will run the "computeMomenta()" function for every model in the tuple. This function might set
 *          calculate density and velocity based on the distributions, for instance.
 */
template<class ...TModel>
inline void Algorithm<TModel...>::calculateMomentaStep() { //...

    if constexpr (sizeof...(TModel) != 0) {

        std::apply([](TModel&... models) {

            (models.computeMomenta(), ...);

        }, mt_Models);

    }

}
