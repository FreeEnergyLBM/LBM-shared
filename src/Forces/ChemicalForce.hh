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

template<class method=Guo, template<class,int> class gradienttype=Gradient>
class ChemicalForce : public ForceBase<method> {
    
    public:

        template<class traits>
        inline double computeXYZ(const int xyz, const int k) const; //Return force at traits::Lattice point k in direction xyz

        template<class traits>
        inline double computeQ(const int xyz, const int k) const;

        template<class traits>
        inline double computeVelocitySource(const int xyz,const int k) const; //Calculate any possible source/correction term for
                                                           //velocity

        template<class traits, int directions>
        inline double computeChemicalForce(int idx, int k) const ;

};

template<class method, template<class,int> class gradienttype>
template<class traits>
inline double ChemicalForce<method, gradienttype>::computeXYZ(const int xyz, const int k) const {
    //tupletype::_;

    if constexpr (has_type<Cartesian,typename method::mt_Stencils>::type::value){
        
        return computeChemicalForce<traits,traits::Lattice::m_NDIM>(xyz,k);
    }
    return 0;

}

template<class method, template<class,int> class gradienttype>
template<class traits>
inline double ChemicalForce<method, gradienttype>::computeQ(const int idx, const int k) const {
    
    if constexpr (has_type<AllDirections,typename method::mt_Stencils>::type::value){
        
        return computeChemicalForce<traits,traits::Stencil::Q>(idx,k);
    }
    return 0;

}

template<class method, template<class,int> class gradienttype>
template<class traits, int directions>
inline double ChemicalForce<method, gradienttype>::computeChemicalForce(int idx, int k) const {
    double sum=0;
    for (int component = 0; component < traits::NumberOfComponents-1; component++) {
        
        sum+=ChemicalPotential<traits::NumberOfComponents-1>::template get<typename traits::Lattice>(k,component) * gradienttype<OrderParameter<traits::NumberOfComponents-1>,(traits::NumberOfComponents-1)>::template get<typename traits::Lattice,directions>(k,component,idx);
        
    }
    return sum;
}

template<class method, template<class,int> class gradienttype>
template<class traits>
inline double ChemicalForce<method, gradienttype>::computeVelocitySource(const int xyz, const int k) const{ //Need to correct velocity
    
    return +computeXYZ<traits>(xyz,k) * traits::Lattice::m_DT / (2.0);
    
}