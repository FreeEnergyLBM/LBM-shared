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

class ChemicalPotentialCalculatorRho : public AddOnBase {
    
    public:

        template<typename traits>
        inline void compute(int k);

        inline void setA(const double A);

        inline void setKappa(const double kappa);

        inline void setRhoLiquid(const double rho) {m_RhoLiquid=rho;};

        inline void setRhoGas(const double rho) {m_RhoVapor=rho;};

        inline void setRho(const double rhol,const double rhov) {m_RhoVapor=rhov;m_RhoLiquid=rhol;};

        double m_A;

        double m_Kappa;

        double m_RhoLiquid;

        double m_RhoVapor;
    
};

template<typename T_traits>
inline void ChemicalPotentialCalculatorRho::compute(int k){

    using Lattice = typename T_traits::Lattice;

    ChemicalPotential<>::get<Lattice>(k) = 2 * m_A * (Density<>::get<Lattice>(k) - m_RhoLiquid)
                                                   * (Density<>::get<Lattice>(k) - m_RhoVapor) 
                                                   * (2*Density<>::get<Lattice>(k) - m_RhoLiquid - m_RhoVapor)
                                           - m_Kappa * LaplacianDensity<>::get<Lattice>(k);

}

inline void ChemicalPotentialCalculatorRho::setA(double A){

    m_A = A;

}

inline void ChemicalPotentialCalculatorRho::setKappa(double kappa){

    m_Kappa = kappa;

}

class ChemicalPotentialCalculatorBinaryLee : public AddOnBase {
    
    public:

        template<typename traits>
        inline void compute(int k);

        inline void setA(double A);

        inline void setKappa(double kappa);

        double m_A;

        double m_Kappa;
    
};

template<typename T_traits>
inline void ChemicalPotentialCalculatorBinaryLee::compute(int k){

    using Lattice = typename T_traits::Lattice;

    ChemicalPotential<>::get<Lattice>(k) = 2 * m_A * OrderParameter<>::get<Lattice>(k)
                                           - 6 * m_A * pow(OrderParameter<>::get<Lattice>(k), 2)
                                           + 4 * m_A * pow(OrderParameter<>::get<Lattice>(k), 3)
                                           - m_Kappa * LaplacianOrderParameter<>::get<Lattice>(k);

}

inline void ChemicalPotentialCalculatorBinaryLee::setA(double A){

    m_A = A;

}

inline void ChemicalPotentialCalculatorBinaryLee::setKappa(double kappa){

    m_Kappa = kappa;

}

class ChemicalPotentialCalculatorNComponent : public AddOnBase {
    
    public:

        double **ma_Gamma;

        double **ma_Beta;

        template<class traits>
        inline void compute(const int k); //Perform any neccessary computations before force is computed

        inline void setA(double **A) {ma_Beta=A;};

        inline void setKappa(double **kappa) {ma_Gamma=kappa;};
    
};

template<class T_traits>
inline void ChemicalPotentialCalculatorNComponent::compute(const int k){ // THIS IS WRONG, NEED 1 - OTHER LAPLACIANS FOR THE FINAL SUM
        
        double gammalaplaciansum=0;
        double sumc=0;

        for (int i=0;i<T_traits::NumberOfComponents-1;i++){
            const double& ci=OrderParameter<T_traits::NumberOfComponents-1>::template get<typename T_traits::Lattice>(k,i);
            double chempot=0;
            sumc=0;
            gammalaplaciansum=0;
            for (int j=0;j<T_traits::NumberOfComponents-1;j++){
                const double& cj=OrderParameter<T_traits::NumberOfComponents-1>::template get<typename T_traits::Lattice>(k,j);
                sumc+=cj;
                const double gammalaplacian = ma_Gamma[i][j]*LaplacianOrderParameter<T_traits::NumberOfComponents-1>::template get<typename T_traits::Lattice,T_traits::NumberOfComponents-1>(k,j);
                gammalaplaciansum += gammalaplacian;
                if (i!=j){
                    
                    chempot+=2*ma_Beta[i][j]*(-12*ci*ci*cj-12*ci*cj*cj+12*ci*cj-4*cj*cj*cj+6*cj*cj-2*cj) - gammalaplacian;
                }
            }
            const double cj=1-sumc;
            
            chempot+=ma_Beta[i][T_traits::NumberOfComponents-1]*(-12*ci*ci*cj-12*ci*cj*cj+12*ci*cj-4*cj*cj*cj+6*cj*cj-2*cj) + gammalaplaciansum;
            ChemicalPotential<T_traits::NumberOfComponents>::template get<typename T_traits::Lattice>(k,i) = chempot;
            
        }
        int i = T_traits::NumberOfComponents-1;
        const double& ci=1-sumc;
        double chempot=0;
        sumc=0;
        gammalaplaciansum=0;
        for (int j=0;j<T_traits::NumberOfComponents-1;j++){
            const double& cj=OrderParameter<T_traits::NumberOfComponents-1>::template get<typename T_traits::Lattice>(k,j);
            sumc+=cj;
            const double gammalaplacian = ma_Gamma[i][j]*LaplacianOrderParameter<T_traits::NumberOfComponents-1>::template get<typename T_traits::Lattice,T_traits::NumberOfComponents-1>(k,j);
            gammalaplaciansum += gammalaplacian;
            if (i!=j){
                
                chempot+=2*ma_Beta[i][j]*(-12*ci*ci*cj-12*ci*cj*cj+12*ci*cj-4*cj*cj*cj+6*cj*cj-2*cj) - gammalaplacian;
            }
        }

        ChemicalPotential<T_traits::NumberOfComponents>::template get<typename T_traits::Lattice>(k,i) = chempot;
        
}