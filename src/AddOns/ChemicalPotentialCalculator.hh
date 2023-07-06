#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "../Service.hh"
#include "AddOnBase.hh"
#include<iostream>
#include<utility>
#include<math.h>

class ChemicalPotentialCalculatorBinary : public AddOnBase {
    
    public:

        template<typename traits>
        inline void compute(int k);

        inline void setA(double A);

        inline void setKappa(double kappa);

        double m_A;

        double m_Kappa;
    
};

template<typename T_traits>
inline void ChemicalPotentialCalculatorBinary::compute(int k){

    using Lattice = typename T_traits::Lattice;

    ChemicalPotential<>::get<Lattice>(k) = -m_A * OrderParameter<>::get<Lattice>(k)
                                                            + m_A * pow(OrderParameter<>::get<Lattice>(k), 3)
                                                            - m_Kappa * LaplacianOrderParameter<>::get<Lattice>(k);

}

inline void ChemicalPotentialCalculatorBinary::setA(double A){

    m_A = A;

}

inline void ChemicalPotentialCalculatorBinary::setKappa(double kappa){

    m_Kappa = kappa;

}