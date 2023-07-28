#pragma once
#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "../Service.hh"
#include "ForceBase.hh"
#include <math.h>
#include<iostream>
#include<utility>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<class method=AllenCahnSourceMethod,int componentID=0>
class AllenCahnSource : public ForceBase<method> {
    
    public:

        template<class traits>
        inline double computeXYZ(int xyz, int k) const; //Return force at traits::Lattice point k in direction xyz

        template<class traits>
        inline double computeQ(int xyz, int k) const;
        
        double m_D;

        double mobility = 0.00333;

        const double magnitudecutoff = 1e-5;

        inline void setAlpha(double alpha){ m_D=alpha; }

        template<class traits>
        inline double computeBeta(int xyz, int k) const;

};

template<class method,int componentID>
template<class traits>
inline double AllenCahnSource<method,componentID>::computeXYZ(int xyz, int k) const {

    double gradx=GradientOrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice,traits::Lattice::NDIM>(k,componentID,0);
    double grady=GradientOrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice,traits::Lattice::NDIM>(k,componentID,1);

    double magnitudegrad2=gradx*gradx+grady*grady;
    if constexpr (traits::Lattice::NDIM==3) magnitudegrad2+=GradientOrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice,traits::Lattice::NDIM>(k,componentID,2)*GradientOrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice,traits::Lattice::NDIM>(k,componentID,2);
    double normal=GradientOrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice,traits::Lattice::NDIM>(k,componentID,xyz)/sqrt(magnitudegrad2);
    
    if (sqrt(magnitudegrad2)>magnitudecutoff) {
        //std::cout<<mobility*(4*OrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice>(k,componentID)*(1.-OrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice>(k,componentID))*normal/m_D);//-computeBeta<traits>(xyz, k))<<std::endl;
        return mobility*(4*OrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice>(k,componentID)*(1.-OrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice>(k,componentID))*normal/m_D-computeBeta<traits>(xyz, k));

    }
    else return 0;
}

template<class method,int componentID>
template<class traits>
inline double AllenCahnSource<method,componentID>::computeQ(int idx, int k) const {

    return 0;

}

template<class method,int componentID>
template<class traits>
inline double AllenCahnSource<method,componentID>::computeBeta(int xyz, int k) const {
    // NORMAL NEEDS TO DEPEND ON COMPONENT
    double sum=0;
    double orderparametersum = 0;
    double gradientsum[traits::Lattice::NDIM]={};
    for (int component = 0; component<traits::NumberOfComponents-1; component++){
            double gradx=GradientOrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice,traits::Lattice::NDIM>(k,component,0);
            gradientsum[0]+=gradx;
            double grady=GradientOrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice,traits::Lattice::NDIM>(k,component,1);
            gradientsum[1]+=grady;
            double magnitudegrad2=gradx*gradx+grady*grady;
            if constexpr (traits::Lattice::NDIM==3) {
                magnitudegrad2+=GradientOrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice,traits::Lattice::NDIM>(k,component,2)*GradientOrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice,traits::Lattice::NDIM>(k,component,2);
                gradientsum[2]+=GradientOrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice,traits::Lattice::NDIM>(k,component,2);
            }
            double normal=GradientOrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice,traits::Lattice::NDIM>(k,component,xyz)/sqrt(magnitudegrad2);
            if (sqrt(magnitudegrad2)>magnitudecutoff) {
                orderparametersum+=OrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice>(k,component);
                sum += 4*OrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice>(k,component)*(1-OrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice>(k,component))*normal/m_D;
            }
    }
    double magnitudegrad2=gradientsum[0]*gradientsum[0]+gradientsum[1]*gradientsum[1];
    if constexpr (traits::Lattice::NDIM==3) magnitudegrad2 += gradientsum[2]*gradientsum[2];
    double normal = (-gradientsum[xyz])/(sqrt(magnitudegrad2));
    if (sqrt(magnitudegrad2)>magnitudecutoff) sum += 4*(1-orderparametersum)*(orderparametersum)*normal/m_D;
    return OrderParameter<traits::NumberOfComponents-1>::template get<typename traits::Lattice>(k,componentID)*sum;

}