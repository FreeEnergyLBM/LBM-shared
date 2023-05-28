#pragma once
#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "AddOnBase.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<class lattice, class method, int num=2>
class ChemicalForceBase : public AddOnBase {
    
    public:

        using Method = method;

        virtual inline double computeXYZ(const int xyz, const int k) const override; //Return force at lattice point k in direction xyz

        virtual inline double computeVelocitySource(const int xyz,const int k) const override; //Calculate any possible source/correction term for
                                                           //velocity

        ChemicalPotential<lattice,num-1> m_ChemicalPotential;

        GradientOrderParameter<lattice,num-1> m_GradOrderParameter;

        LaplacianOrderParameter<lattice,num-1> m_LaplacianOrderParameter;

        OrderParameter<lattice,num-1> m_OrderParameter;

        Density<lattice> m_Density;

};

template<class lattice, class method, int num>
inline double ChemicalForceBase<lattice, method, num>::computeXYZ(const int xyz, const int k) const {

    return m_ChemicalPotential.getParameter(k) * m_GradOrderParameter.getParameterPointer(k)[xyz];

}

template<class lattice, class method, int num>
inline double ChemicalForceBase<lattice, method, num>::computeVelocitySource(const int xyz, const int k) const{ //Need to correct velocity

    return +computeXYZ(xyz,k) * lattice::m_DT / (2.0 * m_Density.getParameter(k));
    
}

template<class lattice, class method>
class ChemicalForceBinary : public ChemicalForceBase<lattice, method, 2> {
    
    public:

        inline void precompute(const int k) override; //Perform any neccessary computations before force is computed

        inline void setA(const double A);

        inline void setKappa(const double kappa);

    private:

        double m_A=0.00015;

        double m_Kappa=0.0003;

};

template<typename lattice, class method>
inline void ChemicalForceBinary<lattice,method>::precompute(const int k){ //Not necessary

    double orderparamcubed=ChemicalForceBase<lattice, method, 2>::m_OrderParameter.getParameter(k) * ChemicalForceBase<lattice, method, 2>::m_OrderParameter.getParameter(k) * ChemicalForceBase<lattice, method, 2>::m_OrderParameter.getParameter(k);
    ChemicalForceBase<lattice, method, 2>::m_ChemicalPotential.getParameter(k) = -m_A * ChemicalForceBase<lattice, method, 2>::m_OrderParameter.getParameter(k) + m_A * orderparamcubed - m_Kappa * ChemicalForceBase<lattice, method, 2>::m_LaplacianOrderParameter.getParameter(k);

}

template<typename lattice, class method>
inline void ChemicalForceBinary<lattice,method>::setA(const double A){ //Not necessary

    m_A=A;

}

template<typename lattice, class method>
inline void ChemicalForceBinary<lattice,method>::setKappa(const double kappa){ //Not necessary

    m_Kappa=kappa;

}