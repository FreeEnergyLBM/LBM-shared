#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "../Service.hh"
#include "AddOnBase.hh"
#include<iostream>
#include<utility>

class ChemicalPotentialCalculatorBinary : public AddOnBase {
    
    public:

        template<typename traits>
        inline void compute(const int k); //Perform any neccessary computations before force is computed

        inline void setA(const double A);

        inline void setKappa(const double kappa);

        double m_A;//=0.00015;

        double m_Kappa;//=0.0003;
    
};

template<typename traits>
inline void ChemicalPotentialCalculatorBinary::compute(const int k){ //Not necessary

    double& orderparam = OrderParameter<>::get<typename traits::Lattice>(k);
    double& chemicalpotential = ChemicalPotential<>::get<typename traits::Lattice>(k);
    double& laplacianorderparameter = LaplacianOrderParameter<>::get<typename traits::Lattice>(k);

    double orderparamcubed =  orderparam * orderparam * orderparam;
    chemicalpotential = -m_A * orderparam + m_A * orderparamcubed - m_Kappa * laplacianorderparameter;

}

inline void ChemicalPotentialCalculatorBinary::setA(const double A){ //Not necessary

    m_A=A;

}

inline void ChemicalPotentialCalculatorBinary::setKappa(const double kappa){ //Not necessary

    m_Kappa=kappa;

}