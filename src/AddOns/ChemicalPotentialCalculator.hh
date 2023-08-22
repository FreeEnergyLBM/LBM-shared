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

        template<class TTraits>
        inline void compute(int k);

        inline void setA(double A);

        inline void setKappa(double kappa);

        double mA;

        double mKappa;
    
};

template<class TTraits>
inline void ChemicalPotentialCalculatorBinary::compute(int k){

    using Lattice = typename TTraits::Lattice;

    ChemicalPotential<>::get<Lattice>(k) = -mA * OrderParameter<>::get<Lattice>(k)
                                                            + mA * pow(OrderParameter<>::get<Lattice>(k), 3)
                                                            - mKappa * LaplacianOrderParameter<>::get<Lattice>(k);

}

inline void ChemicalPotentialCalculatorBinary::setA(double A){

    mA = A;

}

inline void ChemicalPotentialCalculatorBinary::setKappa(double kappa){

    mKappa = kappa;

}

class ChemicalPotentialCalculatorRho : public AddOnBase {
    
    public:

        template<class TTraits>
        inline void compute(int k);

        inline void setA(const double A);

        inline void setKappa(const double kappa);

        inline void setRhoLiquid(const double rho) {mRhoLiquid=rho;};

        inline void setRhoGas(const double rho) {mRhoVapor=rho;};

        inline void setRho(const double rhol,const double rhov) {mRhoVapor=rhov;mRhoLiquid=rhol;};

        double mA;

        double mKappa;

        double mRhoLiquid;

        double mRhoVapor;
    
};

template<class TTraits>
inline void ChemicalPotentialCalculatorRho::compute(int k){

    using Lattice = typename TTraits::Lattice;

    ChemicalPotential<>::get<Lattice>(k) = 2 * mA * (Density<>::get<Lattice>(k) - mRhoLiquid)
                                                   * (Density<>::get<Lattice>(k) - mRhoVapor) 
                                                   * (2*Density<>::get<Lattice>(k) - mRhoLiquid - mRhoVapor)
                                           - mKappa * LaplacianDensity<>::get<Lattice>(k);

}

inline void ChemicalPotentialCalculatorRho::setA(double A){

    mA = A;

}

inline void ChemicalPotentialCalculatorRho::setKappa(double kappa){

    mKappa = kappa;

}

class ChemicalPotentialCalculatorBinaryLee : public AddOnBase {
    
    public:

        template<class TTraits>
        inline void compute(int k);

        inline void setA(double A);

        inline void setKappa(double kappa);

        double mA;

        double mKappa;
    
};

template<class TTraits>
inline void ChemicalPotentialCalculatorBinaryLee::compute(int k){

    using Lattice = typename TTraits::Lattice;

    ChemicalPotential<>::get<Lattice>(k) = 2 * mA * OrderParameter<>::get<Lattice>(k)
                                           - 6 * mA * pow(OrderParameter<>::get<Lattice>(k), 2)
                                           + 4 * mA * pow(OrderParameter<>::get<Lattice>(k), 3)
                                           - mKappa * LaplacianOrderParameter<>::get<Lattice>(k);
    //std::cout<<ChemicalPotential<>::get<Lattice>(k)<<std::endl;
}

inline void ChemicalPotentialCalculatorBinaryLee::setA(double A){

    mA = A;

}

inline void ChemicalPotentialCalculatorBinaryLee::setKappa(double kappa){

    mKappa = kappa;

}

class ChemicalPotentialCalculatorNComponent : public AddOnBase {
    
    public:

        ChemicalPotentialCalculatorNComponent() : mv_Beta({{0}}), mv_Gamma({{0}}) {}

        double **ma_Gamma;

        double **ma_Beta;

        std::vector<std::vector<double>> mv_Beta;
        std::vector<std::vector<double>> mv_Gamma;

        double mD=5;

        template<class traits>
        inline void compute(const int k); //Perform any neccessary computations before force is computed

        inline void setA(double **A) {ma_Beta=A;}
        inline void setKappa(double **kappa) {ma_Gamma=kappa;}
        
        inline void setBeta(int i, int j, double beta) {
            if ((int)mv_Beta.size()-1<i||(int)mv_Beta.size()-1<j){
                mv_Beta.resize(mv_Beta.size()+std::max<int>(i-(int)mv_Beta.size()+1,j-(int)mv_Beta.size()+1));
                for (int l=0;l<(int)mv_Beta.size();l++) {
                    mv_Beta[l].resize(mv_Beta[l].size()+std::max<int>(i-(int)mv_Beta[l].size()+1,j-(int)mv_Beta[l].size()+1));
                    mv_Beta[l][l]=0;
                }
            } 
            if (i!=j) {
                mv_Beta[i][j] = beta;
                mv_Beta[j][i] = beta;
            }            
        }
        inline void setGamma(int i, int j, double gamma) {
            if ((int)mv_Gamma.size()-1<i||(int)mv_Gamma.size()-1<j){
                mv_Gamma.resize(mv_Gamma.size()+std::max<int>(i-(int)mv_Gamma.size()+1,j-(int)mv_Gamma.size()+1));
                for (int l=0;l<(int)mv_Gamma.size();l++) {
                    mv_Gamma[l].resize(mv_Gamma[l].size()+std::max<int>(i-(int)mv_Gamma[l].size()+1,j-(int)mv_Gamma[l].size()+1));
                    mv_Gamma[l][l]=0;
                }
            } 
            if (i!=j) {
                mv_Gamma[i][j] = gamma;
                mv_Gamma[j][i] = gamma;
            }
        }

        inline void setBetaAndGamma(int i, int j, double beta, double gamma) {
            setBeta(i,j,beta);
            setGamma(i,j,gamma);     
        }

        inline void setBeta(std::vector<std::vector<double>>& beta) {
            mv_Beta=beta;
        }
        inline void setGamma(std::vector<std::vector<double>>& gamma) {
            mv_Gamma=gamma;
        }
        inline void setBetaAndGamma(std::vector<std::vector<double>>& beta,std::vector<std::vector<double>>& gamma) {
            mv_Beta=beta;
            mv_Gamma=gamma;
        }

        inline void setD(double d){
            for (int l=0;l<(int)mv_Gamma.size();l++) {
                for (int m=0;m<(int)mv_Gamma.size();m++) {
                    mv_Gamma[l][m]*=d/mD;
                }
            }
            for (int l=0;l<(int)mv_Beta.size();l++) {
                for (int m=0;m<(int)mv_Beta.size();m++) {
                    mv_Beta[l][m]*=mD/d;
                }
            } 
            mD=d;
        }

        inline void setSurfaceTension(int i, int j, double surfacetension) {
            setBeta(i,j,3.0*surfacetension/mD);
            setGamma(i,j,-3.0*mD*surfacetension/4.0);     
        }
    
    public:
        template<int numberofcomponents>
        inline bool checkValid(){
            
            if ((int)mv_Beta.size() != numberofcomponents || (int)mv_Beta[0].size() != numberofcomponents
                || (int)mv_Gamma.size() != numberofcomponents || (int)mv_Gamma[0].size() != numberofcomponents){
                throw std::runtime_error("Number of beta/gamma parameters does not match the number of components.");
                return false;
            }
            return true;
        }

        bool mIsValid;

};

template<class TTraits>
inline void ChemicalPotentialCalculatorNComponent::compute(const int k){ // THIS IS WRONG, NEED - OTHER LAPLACIANS FOR THE FINAL SUM
        
        [[maybe_unused]] static bool isvalid = checkValid<TTraits::NumberOfComponents>();
        mIsValid=isvalid;

        double gammalaplaciansum=0;
        double sumc=0;

        for (int i=0;i<TTraits::NumberOfComponents-1;i++){
            const double& ci=OrderParameter<TTraits::NumberOfComponents-1>::template get<class TTraits::Lattice>(k,i);
            double chempot=0;
            sumc=0;
            gammalaplaciansum=0;
            for (int j=0;j<TTraits::NumberOfComponents-1;j++){
                
                const double& cj=OrderParameter<TTraits::NumberOfComponents-1>::template get<class TTraits::Lattice>(k,j);
                sumc+=cj;
                
                const double gammalaplacian = mv_Gamma[i][j]*LaplacianOrderParameter<TTraits::NumberOfComponents-1>::template get<class TTraits::Lattice,TTraits::NumberOfComponents-1>(k,j);
                
                gammalaplaciansum += gammalaplacian;
                if (i!=j){
                    
                    chempot+=2*mv_Beta[i][j]*(-12*ci*ci*cj-12*ci*cj*cj+12*ci*cj-4*cj*cj*cj+6*cj*cj-2*cj) - gammalaplacian;
                }
                
            }
            const double cj=1-sumc;
            
            chempot+=mv_Beta[i][TTraits::NumberOfComponents-1]*(-12*ci*ci*cj-12*ci*cj*cj+12*ci*cj-4*cj*cj*cj+6*cj*cj-2*cj) + gammalaplaciansum;
            ChemicalPotential<TTraits::NumberOfComponents>::template get<class TTraits::Lattice>(k,i) = chempot;
            
            
        }
        int i = TTraits::NumberOfComponents-1;
        const double& ci=1-sumc;
        double chempot=0;
        sumc=0;
        gammalaplaciansum=0;
        
        for (int j=0;j<TTraits::NumberOfComponents-1;j++){
            const double& cj=OrderParameter<TTraits::NumberOfComponents-1>::template get<class TTraits::Lattice>(k,j);
            sumc+=cj;
            
            const double gammalaplacian = mv_Gamma[i][j]*LaplacianOrderParameter<TTraits::NumberOfComponents-1>::template get<class TTraits::Lattice,TTraits::NumberOfComponents-1>(k,j);
            gammalaplaciansum += gammalaplacian;
            
            if (i!=j){
                
                chempot+=2*mv_Beta[i][j]*(-12*ci*ci*cj-12*ci*cj*cj+12*ci*cj-4*cj*cj*cj+6*cj*cj-2*cj) - gammalaplacian;
            }
        }

        ChemicalPotential<TTraits::NumberOfComponents>::template get<class TTraits::Lattice>(k,i) = chempot;
        
}