#pragma once
#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "../Service.hh"
#include "../Forcing.hh"
#include "ForceBase.hh"
#include<iostream>
#include<utility>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<class T_method=Guo, template<class,int> class T_gradienttype=Gradient>
class ChemicalForceBinary : public ForceBase<T_method> {
    
    public:

        template<class T_traits>
        inline double computeXYZ(int xyz, int k); //Return force at traits::Lattice point k in direction xyz

        template<class T_traits>
        inline double computeQ(int xyz, int k);

        template<class T_traits>
        inline double computeVelocitySource(int xyz, int k); //Calculate any possible source/correction term for velocity

        template<class T_traits, int T_directions>
        inline double computeChemicalForce(int idx, int k);

};

template<class T_method, template<class,int> class T_gradienttype>
template<class T_traits>
inline double ChemicalForceBinary<T_method, T_gradienttype>::computeXYZ(int xyz, int k) {
    //tupletype::_;

    if constexpr (has_type<Cartesian,typename T_method::mt_Stencils>::type::value){
        
        return computeChemicalForce<T_traits,T_traits::Lattice::NDIM>(xyz,k);

    }
    return 0;

}

template<class T_method, template<class,int> class T_gradienttype>
template<class T_traits>
inline double ChemicalForceBinary<T_method, T_gradienttype>::computeQ(int idx, int k) {
    
    if constexpr (has_type<AllDirections,typename T_method::mt_Stencils>::type::value){
        
        return computeChemicalForce<T_traits,T_traits::Stencil::Q>(idx,k);

    }
    return 0;

}

template<class T_method, template<class,int> class T_gradienttype>
template<class T_traits, int T_directions>
inline double ChemicalForceBinary<T_method, T_gradienttype>::computeChemicalForce(int idx, int k) {

    double sum = 0;

    for (int component = 0; component < T_traits::NumberOfComponents - 1; component++) {
        
        sum += ChemicalPotential<T_traits::NumberOfComponents - 1>::template get<typename T_traits::Lattice>(k, component) * T_gradienttype<OrderParameter<T_traits::NumberOfComponents - 1>,(T_traits::NumberOfComponents - 1)>::template get<typename T_traits::Lattice, T_directions>(k, component, idx);
        
    }

    return sum;

}

template<class T_method=Guo, template<class,int> class T_gradienttype=Gradient>
class ChemicalForce : public ForceBase<T_method> {
    
    public:

        template<class T_traits>
        inline double computeXYZ(int xyz, int k); //Return force at traits::Lattice point k in direction xyz

        template<class T_traits>
        inline double computeQ(int xyz, int k);

        template<class T_traits>
        inline double computeVelocitySource(int xyz, int k); //Calculate any possible source/correction term for velocity

        template<class T_traits, int T_directions>
        inline double computeChemicalForce(int idx, int k);

};

template<class T_method, template<class,int> class T_gradienttype>
template<class T_traits>
inline double ChemicalForce<T_method, T_gradienttype>::computeXYZ(int xyz, int k) {
    //tupletype::_;

    if constexpr (has_type<Cartesian,typename T_method::mt_Stencils>::type::value){
        
        return computeChemicalForce<T_traits,T_traits::Lattice::NDIM>(xyz,k);

    }
    return 0;

}

template<class T_method, template<class,int> class T_gradienttype>
template<class T_traits>
inline double ChemicalForce<T_method, T_gradienttype>::computeQ(int idx, int k) {
    
    if constexpr (has_type<AllDirections,typename T_method::mt_Stencils>::type::value){
        
        return computeChemicalForce<T_traits,T_traits::Stencil::Q>(idx,k);

    }
    return 0;

}

template<class T_method, template<class,int> class T_gradienttype>
template<class T_traits, int T_directions>
inline double ChemicalForce<T_method, T_gradienttype>::computeChemicalForce(int idx, int k) {

    double sum = 0;

    for (int component = 0; component < T_traits::NumberOfComponents - 1; component++) {
        
        sum += ChemicalPotential<T_traits::NumberOfComponents>::template get<typename T_traits::Lattice>(k, component) * T_gradienttype<OrderParameter<T_traits::NumberOfComponents - 1>,(T_traits::NumberOfComponents - 1)>::template get<typename T_traits::Lattice, T_directions>(k, component, idx);
        
    }

    return sum;

}

template<class T_method, template<class,int> class T_gradienttype>
template<class T_traits>
inline double ChemicalForce<T_method, T_gradienttype>::computeVelocitySource(const int xyz, const int k) { //Need to correct velocity
    
    return +computeXYZ<T_traits>(xyz, k) * T_traits::Lattice::DT / (2.0);
    
}


template<class T_method=Guo, template<class,int> class T_gradienttype=Gradient>
class ChemicalForceRho : public ForceBase<T_method> {
    
    public:

        template<class T_traits>
        inline double computeXYZ(int xyz, int k); //Return force at traits::Lattice point k in direction xyz

        template<class T_traits>
        inline double computeQ(int xyz, int k);

        template<class T_traits>
        inline double computeVelocitySource(int xyz, int k); //Calculate any possible source/correction term for velocity

        template<class T_traits, int T_directions>
        inline double computeChemicalForce(int idx, int k);

};

template<class T_method, template<class,int> class T_gradienttype>
template<class T_traits>
inline double ChemicalForceRho<T_method, T_gradienttype>::computeXYZ(int xyz, int k) {
    //tupletype::_;

    if constexpr (has_type<Cartesian,typename T_method::mt_Stencils>::type::value){
        
        return computeChemicalForce<T_traits,T_traits::Lattice::NDIM>(xyz,k);

    }
    return 0;

}

template<class T_method, template<class,int> class T_gradienttype>
template<class T_traits>
inline double ChemicalForceRho<T_method, T_gradienttype>::computeQ(int idx, int k) {
    
    if constexpr (has_type<AllDirections,typename T_method::mt_Stencils>::type::value){
        
        return computeChemicalForce<T_traits,T_traits::Stencil::Q>(idx,k);

    }
    return 0;

}

template<class T_method, template<class,int> class T_gradienttype>
template<class T_traits, int T_directions>
inline double ChemicalForceRho<T_method, T_gradienttype>::computeChemicalForce(int idx, int k) {

    double sum = 0;

    for (int component = 0; component < T_traits::NumberOfComponents - 1; component++) {
        
        sum += ChemicalPotential<T_traits::NumberOfComponents - 1>::template get<typename T_traits::Lattice>(k, component) * T_gradienttype<Density<T_traits::NumberOfComponents - 1>,(T_traits::NumberOfComponents - 1)>::template get<typename T_traits::Lattice, T_directions>(k, component, idx);
        
    }

    return sum;

}

template<class T_method, template<class,int> class T_gradienttype>
template<class T_traits>
inline double ChemicalForceRho<T_method, T_gradienttype>::computeVelocitySource(const int xyz, const int k) { //Need to correct velocity
    
    return +computeXYZ<T_traits>(xyz, k) * T_traits::Lattice::DT / (2.0);
    
}