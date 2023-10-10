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

template<class TMethod=Guo, template<class,int> class TGradientType=Gradient>
class ChemicalForceBinary : public ForceBase<TMethod> {
    
    public:

        template<class TTraits>
        inline double computeXYZ(int xyz, int k); //Return force at traits::Lattice point k in direction xyz

        template<class TTraits>
        inline double computeQ(int xyz, int k);

        template<class TTraits>
        inline double computeVelocitySource(int xyz, int k); //Calculate any possible source/correction term for velocity

        template<class TTraits, int TDirections>
        inline double computeChemicalForce(int idx, int k);

};

template<int i>
void test(){}

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits>
inline double ChemicalForceBinary<TMethod, TGradientType>::computeXYZ(int xyz, int k) {
    //tupletype::_;

    if constexpr (has_type<Cartesian,typename TMethod::mt_Stencils>::type::value){
        
        return computeChemicalForce<TTraits,TTraits::Lattice::NDIM>(xyz,k);
        
    }
    return 0;

}

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits>
inline double ChemicalForceBinary<TMethod, TGradientType>::computeQ(int idx, int k) {
    
    if constexpr (has_type<AllDirections,typename TMethod::mt_Stencils>::type::value){
        
        return computeChemicalForce<TTraits,TTraits::Stencil::Q>(idx,k);

    }
    return 0;

}

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits, int TDirections>
inline double ChemicalForceBinary<TMethod, TGradientType>::computeChemicalForce(int idx, int k) {

    double sum = 0;

    for (int component = 0; component < TTraits::NumberOfComponents - 1; component++) {
        //std::cout<<ChemicalPotential<>::template get<typename TTraits::Lattice>(k)<<std::endl;
        double chemPot = ChemicalPotential<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice>(k, component);
        double gradOP = TGradientType<OrderParameter<TTraits::NumberOfComponents - 1>,(TTraits::NumberOfComponents - 1)>::template get<typename TTraits::Lattice, TDirections>(k, component, idx);
        sum += chemPot * gradOP;
        //
    }
    
    return sum;

}

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits>
inline double ChemicalForceBinary<TMethod, TGradientType>::computeVelocitySource(const int xyz, const int k) { //Need to correct velocity
    
    return +computeXYZ<TTraits>(xyz, k) * TTraits::Lattice::DT / (2.0);
    
}


template<class TMethod=Guo, template<class,int> class TGradientType=Gradient>
class ChemicalForce : public ForceBase<TMethod> {
    
    public:

        template<class TTraits>
        inline double computeXYZ(int xyz, int k); //Return force at traits::Lattice point k in direction xyz

        template<class TTraits>
        inline double computeQ(int xyz, int k);

        template<class TTraits>
        inline double computeVelocitySource(int xyz, int k); //Calculate any possible source/correction term for velocity

        template<class TTraits, int TDirections>
        inline double computeChemicalForce(int idx, int k);

};

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits>
inline double ChemicalForce<TMethod, TGradientType>::computeXYZ(int xyz, int k) {
    //tupletype::_;

    if constexpr (has_type<Cartesian,typename TMethod::mt_Stencils>::type::value){
        
        return computeChemicalForce<TTraits,TTraits::Lattice::NDIM>(xyz,k);

    }
    return 0;

}

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits>
inline double ChemicalForce<TMethod, TGradientType>::computeQ(int idx, int k) {
    
    if constexpr (has_type<AllDirections,typename TMethod::mt_Stencils>::type::value){
        
        return computeChemicalForce<TTraits,TTraits::Stencil::Q>(idx,k);

    }
    return 0;

}

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits, int TDirections>
inline double ChemicalForce<TMethod, TGradientType>::computeChemicalForce(int idx, int k) {

    double sum = 0;

    double gradopsum = 0;

    for (int component = 0; component < TTraits::NumberOfComponents - 1; component++) {
        
        const double& chemPot = ChemicalPotential<TTraits::NumberOfComponents>::template get<typename TTraits::Lattice>(k, component);
        const double& gradOP = TGradientType<OrderParameter<TTraits::NumberOfComponents - 1>,(TTraits::NumberOfComponents - 1)>::template get<typename TTraits::Lattice, TDirections>(k, component, idx);
        sum += chemPot * gradOP;
        gradopsum += gradOP;

    }

    sum += ChemicalPotential<TTraits::NumberOfComponents>::template get<typename TTraits::Lattice>(k, TTraits::NumberOfComponents - 1) * (-gradopsum);

    return sum;

}

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits>
inline double ChemicalForce<TMethod, TGradientType>::computeVelocitySource(const int xyz, const int k) { //Need to correct velocity
    
    return +computeXYZ<TTraits>(xyz, k) * TTraits::Lattice::DT / (2.0);
    
}


template<class TMethod=Guo, template<class,int> class TGradientType=Gradient>
class ChemicalForceRho : public ForceBase<TMethod> {
    
    public:

        template<class TTraits>
        inline double computeXYZ(int xyz, int k); //Return force at traits::Lattice point k in direction xyz

        template<class TTraits>
        inline double computeQ(int xyz, int k);

        template<class TTraits>
        inline double computeVelocitySource(int xyz, int k); //Calculate any possible source/correction term for velocity

        template<class TTraits, int TDirections>
        inline double computeChemicalForce(int idx, int k);

};

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits>
inline double ChemicalForceRho<TMethod, TGradientType>::computeXYZ(int xyz, int k) {
    //tupletype::_;

    if constexpr (has_type<Cartesian,typename TMethod::mt_Stencils>::type::value){
        
        return computeChemicalForce<TTraits,TTraits::Lattice::NDIM>(xyz,k);

    }
    return 0;

}

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits>
inline double ChemicalForceRho<TMethod, TGradientType>::computeQ(int idx, int k) {
    
    if constexpr (has_type<AllDirections,typename TMethod::mt_Stencils>::type::value){
        
        return computeChemicalForce<TTraits,TTraits::Stencil::Q>(idx,k);

    }
    return 0;

}

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits, int TDirections>
inline double ChemicalForceRho<TMethod, TGradientType>::computeChemicalForce(int idx, int k) {

    double sum = 0;

    for (int component = 0; component < TTraits::NumberOfComponents - 1; component++) {
        
        sum += ChemicalPotential<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice>(k, component) * TGradientType<Density<TTraits::NumberOfComponents - 1>,(TTraits::NumberOfComponents - 1)>::template get<typename TTraits::Lattice, TDirections>(k, component, idx);
        
    }

    return sum;

}

template<class TMethod, template<class,int> class TGradientType>
template<class TTraits>
inline double ChemicalForceRho<TMethod, TGradientType>::computeVelocitySource(const int xyz, const int k) { //Need to correct velocity
    
    return +computeXYZ<TTraits>(xyz, k) * TTraits::Lattice::DT / (2.0);
    
}
