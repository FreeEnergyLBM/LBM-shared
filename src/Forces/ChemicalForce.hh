#pragma once
#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "ForceBase.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<class lattice>
class ChemicalForce : public ForceBase {
    
    public:

        inline double compute(const int xyz, const int k) const override; //Return force at lattice point k in direction xyz

        inline void precompute(const int k) override; //Perform any neccessary computations before force is computed

        inline double computeVelocitySource(const int xyz,const int k) const override; //Calculate any possible source/correction term for
                                                           //velocity

    private:

        double m_A=0.00015;

        double m_Kappa=0.0003;

        ChemicalPotential<lattice> m_ChemicalPotential;

        GradientOrderParameter<lattice> m_GradOrderParameter;

        LaplacianOrderParameter<lattice> m_LaplacianOrderParameter;

        OrderParameter<lattice> m_OrderParameter;

        Density<lattice> m_Density;

};

template<typename lattice>
inline double ChemicalForce<lattice>::compute(const int xyz, const int k) const {

    return m_ChemicalPotential.getParameter(k) * m_GradOrderParameter.getParameterPointer(k)[xyz];

}

template<typename lattice>
inline void ChemicalForce<lattice>::precompute(const int k){ //Not necessary

    double orderparamcubed=m_OrderParameter.getParameter(k) * m_OrderParameter.getParameter(k) * m_OrderParameter.getParameter(k);
    m_ChemicalPotential.getParameter(k) = -m_A * m_OrderParameter.getParameter(k) + m_A * orderparamcubed - m_Kappa * m_LaplacianOrderParameter.getParameter(k);

}

template<typename lattice>
inline double ChemicalForce<lattice>::computeVelocitySource(const int xyz, const int k) const{ //Need to correct velocity

    return +compute(xyz,k) * lattice::m_DT / (2.0 * m_Density.getParameter(k));
    
}