#pragma once
#include "../Geometry.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "../Service.hh"
#include "AddOnBase.hh"

class ChemicalPotentialCalculatorBinary : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

    inline void setA(double A);

    inline void setKappa(double kappa);

    double mA;
    double mKappa;
};

template <class TTraits>
inline void ChemicalPotentialCalculatorBinary::compute(int k) {
    using Lattice = typename TTraits::Lattice;

    if (Geometry<Lattice>::isBulkSolid(k)) return;

    const double& orderparameter = OrderParameter<>::get<Lattice>(k);
    const double& laplacian = LaplacianOrderParameter<>::get<Lattice>(k);

    ChemicalPotential<>::get<Lattice>(k) = -mA * orderparameter + mA * pow(orderparameter, 3) - mKappa * laplacian;
}

inline void ChemicalPotentialCalculatorBinary::setA(double A) { mA = A; }

inline void ChemicalPotentialCalculatorBinary::setKappa(double kappa) { mKappa = kappa; }

class ChemicalPotentialCalculatorRho : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

    inline void setA(const double A);

    inline void setKappa(const double kappa);

    inline void setRhoLiquid(const double rho) { mRhoLiquid = rho; };

    inline void setRhoGas(const double rho) { mRhoVapor = rho; };

    inline void setRho(const double rhol, const double rhov) {
        mRhoVapor = rhov;
        mRhoLiquid = rhol;
    };

    double mA;
    double mKappa;
    double mRhoLiquid;
    double mRhoVapor;
};

template <class TTraits>
inline void ChemicalPotentialCalculatorRho::compute(int k) {
    using Lattice = typename TTraits::Lattice;

    if (Geometry<Lattice>::isBulkSolid(k)) return;

    const double& density = Density<>::get<Lattice>(k);
    const double& laplacian = LaplacianDensity<>::get<Lattice>(k);

    ChemicalPotential<>::get<Lattice>(k) =
        2 * mA * (density - mRhoLiquid) * (density - mRhoVapor) * (2 * density - mRhoLiquid - mRhoVapor) -
        mKappa * laplacian;
}

inline void ChemicalPotentialCalculatorRho::setA(double A) { mA = A; }

inline void ChemicalPotentialCalculatorRho::setKappa(double kappa) { mKappa = kappa; }

class ChemicalPotentialCalculatorBinaryLee : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

    template <class TTraits>
    inline void communicate() {
        using Lattice = typename TTraits::Lattice;
        Lattice::communicate(ChemicalPotential<>::getInstance<Lattice>());
    }

    inline void setA(double A);

    inline void setKappa(double kappa);

    inline void setOmega(double omega);

    double mA;
    double mKappa;
    double mOmega = 0.0001;
};

template <class TTraits>
inline void ChemicalPotentialCalculatorBinaryLee::compute(int k) {
    using Lattice = typename TTraits::Lattice;

    if (Geometry<Lattice>::isBulkSolid(k)) return;

    const double& orderparameter = OrderParameter<>::get<Lattice>(k);
    const double& laplacian = LaplacianOrderParameter<>::get<Lattice>(k);

    ChemicalPotential<>::get<Lattice>(k) = 2 * mA * orderparameter - 6 * mA * pow(orderparameter, 2) +
                                           4 * mA * pow(orderparameter, 3) +
                                           2 * mOmega * (orderparameter < 0) * orderparameter - mKappa * laplacian;
}

inline void ChemicalPotentialCalculatorBinaryLee::setA(double A) { mA = A; }

inline void ChemicalPotentialCalculatorBinaryLee::setKappa(double kappa) { mKappa = kappa; }

inline void ChemicalPotentialCalculatorBinaryLee::setOmega(double omega) { mOmega = omega; }

class ChemicalPotentialCalculatorTernaryLee : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

    inline void setSurfaceTension(double sigma12, double sigma13, double sigma23);

    inline void setOmega(double omega);

    inline void setLambda(double lambda);

    inline void setInterfaceWidth(double interfacewidth);

    template <class TTraits>
    inline void communicate() {
        using Lattice = typename TTraits::Lattice;
        Lattice::communicate(ChemicalPotential<0>::getInstance<Lattice>());
        Lattice::communicate(ChemicalPotential<1>::getInstance<Lattice>());
        Lattice::communicate(ChemicalPotential<2>::getInstance<Lattice>());
    }

    std::array<double, 3> mSigma;
    std::array<double, 3> mGamma;
    double mGammaT;
    double mOmega = 0.0001;
    double mLambda = 0;
    double mInterfaceWidth;
};

template <class TTraits>
inline void ChemicalPotentialCalculatorTernaryLee::compute(int k) {
    using Lattice = typename TTraits::Lattice;

    if (Geometry<Lattice>::isBulkSolid(k)) return;

    const double C1 = OrderParameter<0>::get<Lattice>(k);  // C1, C2, C3
    const double C2 = OrderParameter<1>::get<Lattice>(k);
    const double C3 = 1 - C1 - C2;
    using Lap1 = LaplacianOrderParameter<0>;
    using Lap2 = LaplacianOrderParameter<1>;

    double dEd12 = mGamma[0] * C1 * (1 - 3 * C1 + 2 * C1 * C1) - mGamma[1] * C2 * (1 - 3 * C2 + 2 * C2 * C2) +
                   2 * mLambda * C3 * C3 * (C2 * C2 * C1 - C1 * C1 * C2) + (C1 < 0) * 2 * mOmega * C1 -
                   (C2 < 0) * 2 * mOmega * C2;

    double dEd13 = mGamma[0] * C1 * (1 - 3 * C1 + 2 * C1 * C1) - mGamma[2] * C3 * (1 - 3 * C3 + 2 * C3 * C3) +
                   2 * mLambda * C2 * C2 * (C3 * C3 * C1 - C1 * C1 * C3) + (C1 < 0) * 2 * mOmega * C1 -
                   (C3 < 0) * 2 * mOmega * C3;

    double dEd23 = mGamma[1] * C2 * (1 - 3 * C2 + 2 * C2 * C2) - mGamma[2] * C3 * (1 - 3 * C3 + 2 * C3 * C3) +
                   2 * mLambda * C1 * C1 * (C3 * C3 * C2 - C2 * C2 * C3) + (C2 < 0) * 2 * mOmega * C2 -
                   (C3 < 0) * 2 * mOmega * C3;

    ChemicalPotential<0>::get<Lattice>(k) =
        (4.0 * mGammaT / mInterfaceWidth * (1 / mGamma[1] * dEd12 + 1 / mGamma[2] * dEd13) -
         3.0 / 4.0 * mInterfaceWidth * mGamma[0] * (Lap1::get<Lattice>(k)));

    ChemicalPotential<1>::get<Lattice>(k) =
        (4.0 * mGammaT / mInterfaceWidth * (-1 / mGamma[0] * dEd12 + 1 / mGamma[2] * dEd23) -
         3.0 / 4.0 * mInterfaceWidth * mGamma[1] * (Lap2::get<Lattice>(k)));

    ChemicalPotential<2>::get<Lattice>(k) =
        (4.0 * mGammaT / mInterfaceWidth * (-1 / mGamma[0] * dEd13 - 1 / mGamma[1] * dEd23) -
         3.0 / 4.0 * mInterfaceWidth * mGamma[2] * (-Lap1::get<Lattice>(k) - Lap2::get<Lattice>(k)));
}

inline void ChemicalPotentialCalculatorTernaryLee::setSurfaceTension(double sigma12, double sigma13, double sigma23) {
    mSigma[0] = sigma12;
    mSigma[1] = sigma13;
    mSigma[2] = sigma23;

    mGamma[0] = sigma12 + sigma13 - sigma23;
    mGamma[1] = sigma12 + sigma23 - sigma13;
    mGamma[2] = sigma13 + sigma23 - sigma12;

    mGammaT = 3.0 / (1.0 / mGamma[0] + 1.0 / mGamma[1] + 1.0 / mGamma[2]);
}

inline void ChemicalPotentialCalculatorTernaryLee::setOmega(double omega) { mOmega = omega; }

inline void ChemicalPotentialCalculatorTernaryLee::setLambda(double lambda) { mLambda = lambda; }

inline void ChemicalPotentialCalculatorTernaryLee::setInterfaceWidth(double interfacewidth) {
    mInterfaceWidth = interfacewidth;
}

class ChemicalPotentialCalculatorTernaryLeeExtraPotential : public ChemicalPotentialCalculatorTernaryLee {
   public:
    template <class TTraits>
    inline void compute(int k);
    inline void setPreOmega(double preomega) { mPreOmega = preomega; }
    inline void setPotentialComponent(double component) { mComponent = component; }
    double mPreOmega = 0.00333;  // 1;
    inline void setPotentialCondition(double (*condition)(int k)) { evalPotentialCondition = condition; }

   private:
    static double defaultCondition(int k) { return true; }
    double (*evalPotentialCondition)(int k) = &defaultCondition;
    int mComponent = 0;
};

template <class TTraits>
inline void ChemicalPotentialCalculatorTernaryLeeExtraPotential::compute(int k) {
    using Lattice = typename TTraits::Lattice;

    if (Geometry<Lattice>::isBulkSolid(k)) return;

    // std::cout<<ChemicalPotential<>::get<Lattice>(k)<<std::endl;
    const double C1 = OrderParameter<0>::get<Lattice>(k);  // C1, C2, C3
    const double C2 = OrderParameter<1>::get<Lattice>(k);
    const double C3 = 1 - C1 - C2;
    double chem1 = 0.0;
    double chem2 = 0.0;
    double chem3 = 0.0;

    double dEd12 = mGamma[0] * C1 * (1 - 3 * C1 + 2 * C1 * C1) - mGamma[1] * C2 * (1 - 3 * C2 + 2 * C2 * C2) +
                   2 * mLambda * C3 * C3 * (C2 * C2 * C1 - C1 * C1 * C2) + (C1 < 0) * 2 * mOmega * C1 -
                   (C2 < 0) * 2 * mOmega * C2;

    double dEd13 = mGamma[0] * C1 * (1 - 3 * C1 + 2 * C1 * C1) - mGamma[2] * C3 * (1 - 3 * C3 + 2 * C3 * C3) +
                   2 * mLambda * C2 * C2 * (C3 * C3 * C1 - C1 * C1 * C3) + (C1 < 0) * 2 * mOmega * C1 -
                   (C3 < 0) * 2 * mOmega * C3;

    double dEd23 = mGamma[1] * C2 * (1 - 3 * C2 + 2 * C2 * C2) - mGamma[2] * C3 * (1 - 3 * C3 + 2 * C3 * C3) +
                   2 * mLambda * C1 * C1 * (C3 * C3 * C2 - C2 * C2 * C3) + (C2 < 0) * 2 * mOmega * C2 -
                   (C3 < 0) * 2 * mOmega * C3;

    //* Extra Potential
    double pre_omega = mPreOmega;

    double conf_pot[3] = {0.0, 0.0, 0.0};
    double pot_sign[3] = {0.0, 0.0, 0.0};
    double Chem_P_extra[3] = {0.0, 0.0, 0.0};
    /*if (evalPotentialCondition(k)) {
        if(OrderParameter<1>::get<Lattice>(k)<0.5){
            conf_pot[mComponent] = -1.0;
            pot_sign[mComponent] = 0.0;
        }
    }
    else if (OrderParameter<1>::get<Lattice>(k)>=0.5) {
        conf_pot[mComponent] = 1.0;
        pot_sign[mComponent] = 0.0;
    }*/ //conf_pot[mComponent] = evalPotentialCondition(k);

    /*Chem_P_extra[0] = pre_omega * 12.0 * conf_pot[0] * (2.0 * pot_sign[0] - 1.0) * (pot_sign[0] - C1) *
                      (pot_sign[0] - C1) * (pot_sign[0] + C1 - 1.0);
    Chem_P_extra[1] = pre_omega * 12.0 * conf_pot[1] * (2.0 * pot_sign[1] - 1.0) * (pot_sign[1] - C2) *
                      (pot_sign[1] - C2) * (pot_sign[1] + C2 - 1.0);
    Chem_P_extra[2] = -pre_omega * 12.0 * conf_pot[1] * (2.0 * pot_sign[2] - 1.0) * (pot_sign[2] - C3) *
                      (pot_sign[2] - C3) * (pot_sign[2] + C3 - 1.0);*/

    // Chem_P_extra[0] = evalPotentialCondition(k) * pre_omega * 6.0 * (C1) *
    //                   (C1) * (2*C1 - 1.0);
    // Chem_P_extra[1] = evalPotentialCondition(k) * pre_omega * 6.0 * (C2) *
    //                   (C2) * (2*C2 - 1.0);
    Chem_P_extra[1] = evalPotentialCondition(k) * pre_omega * (2.0 * C2 - 1.0);
    // Chem_P_extra[2] = evalPotentialCondition(k) * pre_omega * 6.0 * (C3) *
    //                   (C3) * (2*C3 - 1.0);

    dEd12 += Chem_P_extra[0] - Chem_P_extra[1];
    dEd13 += Chem_P_extra[0] - Chem_P_extra[2];
    dEd23 += Chem_P_extra[1] - Chem_P_extra[2];
    // Extra potential

    /*ChemicalPotential2<0>::get<Lattice>(k) = 4.0 * mGammaT / mInterfaceWidth * (1 / mGamma[1] * (Chem_P_extra[0] -
    Chem_P_extra[1]) + 1 / mGamma[2] * (Chem_P_extra[0] - Chem_P_extra[2])); ChemicalPotential2<1>::get<Lattice>(k)
    = 4.0 * mGammaT / mInterfaceWidth * (-1 / mGamma[0] * (Chem_P_extra[0] - Chem_P_extra[1]) + 1 / mGamma[2] *
    (Chem_P_extra[1] - Chem_P_extra[2])); ChemicalPotential2<2>::get<Lattice>(k) = 4.0 * mGammaT / mInterfaceWidth *
    (-1 / mGamma[0] * (Chem_P_extra[0] - Chem_P_extra[2]) - 1 / mGamma[1] * (Chem_P_extra[1] - Chem_P_extra[2]));*/

    chem1 = (4.0 * mGammaT / mInterfaceWidth * (1 / mGamma[1] * dEd12 + 1 / mGamma[2] * dEd13) -
             3.0 / 4.0 * mInterfaceWidth * mGamma[0] * (LaplacianOrderParameter<0>::get<Lattice>(k)));

    chem2 = (4.0 * mGammaT / mInterfaceWidth * (-1 / mGamma[0] * dEd12 + 1 / mGamma[2] * dEd23) -
             3.0 / 4.0 * mInterfaceWidth * mGamma[1] * (LaplacianOrderParameter<1>::get<Lattice>(k)));

    chem3 = (4.0 * mGammaT / mInterfaceWidth * (-1 / mGamma[0] * dEd13 - 1 / mGamma[1] * dEd23) -
             3.0 / 4.0 * mInterfaceWidth * mGamma[2] *
                 (-LaplacianOrderParameter<0>::get<Lattice>(k) - LaplacianOrderParameter<1>::get<Lattice>(k)));

    ChemicalPotential<0>::get<Lattice>(k) = chem1;
    ChemicalPotential<1>::get<Lattice>(k) = chem2;
    ChemicalPotential<2>::get<Lattice>(k) = chem3;
}

class ChemicalPotentialCalculatorNComponent : public AddOnBase {
   public:
    ChemicalPotentialCalculatorNComponent() : mv_Beta({{0}}), mv_Gamma({{0}}) {}

    template <class traits>
    inline void compute(const int k);

    inline void setA(double** A) { ma_Beta = A; }

    inline void setKappa(double** kappa) { ma_Gamma = kappa; }

    inline void setBeta(int i, int j, double beta) {
        if ((int)mv_Beta.size() - 1 < i || (int)mv_Beta.size() - 1 < j) {
            mv_Beta.resize(mv_Beta.size() + std::max<int>(i - (int)mv_Beta.size() + 1, j - (int)mv_Beta.size() + 1));
            for (int l = 0; l < (int)mv_Beta.size(); l++) {
                mv_Beta[l].resize(mv_Beta[l].size() +
                                  std::max<int>(i - (int)mv_Beta[l].size() + 1, j - (int)mv_Beta[l].size() + 1));
                mv_Beta[l][l] = 0;
            }
        }
        if (i != j) {
            mv_Beta[i][j] = beta;
            mv_Beta[j][i] = beta;
        }
    }
    inline void setGamma(int i, int j, double gamma) {
        if ((int)mv_Gamma.size() - 1 < i || (int)mv_Gamma.size() - 1 < j) {
            mv_Gamma.resize(mv_Gamma.size() +
                            std::max<int>(i - (int)mv_Gamma.size() + 1, j - (int)mv_Gamma.size() + 1));
            for (int l = 0; l < (int)mv_Gamma.size(); l++) {
                mv_Gamma[l].resize(mv_Gamma[l].size() +
                                   std::max<int>(i - (int)mv_Gamma[l].size() + 1, j - (int)mv_Gamma[l].size() + 1));
                mv_Gamma[l][l] = 0;
            }
        }
        if (i != j) {
            mv_Gamma[i][j] = gamma;
            mv_Gamma[j][i] = gamma;
        }
    }

    inline void setBetaAndGamma(int i, int j, double beta, double gamma) {
        setBeta(i, j, beta);
        setGamma(i, j, gamma);
    }

    inline void setBeta(std::vector<std::vector<double>>& beta) { mv_Beta = beta; }
    inline void setGamma(std::vector<std::vector<double>>& gamma) { mv_Gamma = gamma; }
    inline void setBetaAndGamma(std::vector<std::vector<double>>& beta, std::vector<std::vector<double>>& gamma) {
        mv_Beta = beta;
        mv_Gamma = gamma;
    }

    inline void setD(double d) {
        for (int l = 0; l < (int)mv_Gamma.size(); l++) {
            for (int m = 0; m < (int)mv_Gamma.size(); m++) {
                mv_Gamma[l][m] *= d / mD;
            }
        }
        for (int l = 0; l < (int)mv_Beta.size(); l++) {
            for (int m = 0; m < (int)mv_Beta.size(); m++) {
                mv_Beta[l][m] *= mD / d;
            }
        }
        mD = d;
    }

    inline void setSurfaceTension(int i, int j, double surfacetension) {
        setBeta(i, j, 3.0 * surfacetension / mD);
        setGamma(i, j, -3.0 * mD * surfacetension / 4.0);
    }

   public:
    template <int numberofcomponents>
    inline bool checkValid() {
        if ((int)mv_Beta.size() != numberofcomponents || (int)mv_Beta[0].size() != numberofcomponents ||
            (int)mv_Gamma.size() != numberofcomponents || (int)mv_Gamma[0].size() != numberofcomponents) {
            throw std::runtime_error(
                "Number of beta/gamma parameters does not match the number of "
                "components.");
            return false;
        }
        return true;
    }

    bool mIsValid;
    double** ma_Gamma;
    double** ma_Beta;
    std::vector<std::vector<double>> mv_Beta;
    std::vector<std::vector<double>> mv_Gamma;
    double mD = 5;
};

template <class TTraits>
inline void ChemicalPotentialCalculatorNComponent::compute(const int k) {
    if (Geometry<typename TTraits::Lattice>::isBulkSolid(k)) return;

    [[maybe_unused]] static bool isvalid = checkValid<TTraits::NumberOfComponents>();
    mIsValid = isvalid;

    using Lattice = typename TTraits::Lattice;
    constexpr int N = TTraits::NumberOfComponents;

    double chempot = 0.0;
    double chempot_bulk = 0.0;
    double chempot_tension = 0.0;
    double sumci = 0.0;
    double sumcj = 0.0;
    double onlylaplacian = 0.0;
    double gammalaplacian = 0.0;
    double onlylaplaciansum = 0.0;
    double gammalaplaciansum = 0.0;

    for (int i = 0; i < N - 1; i++) {
        chempot = 0.0;
        chempot_bulk = 0.0;
        chempot_tension = 0.0;
        const double& ci = getInstance<OrderParameter, N, Lattice>(i)[k];
        sumci += ci;

        sumcj = 0.0;
        gammalaplaciansum = 0.0;
        onlylaplaciansum = 0.0;

        for (int j = 0; j < N - 1; j++) {
            const double& cj = getInstance<OrderParameter, N, Lattice>(j)[k];
            sumcj += cj;

            onlylaplacian = getInstance<LaplacianOrderParameter, N, Lattice>(j)[k];

            onlylaplaciansum += onlylaplacian;
            gammalaplacian = mv_Gamma[i][j] * onlylaplacian;
            gammalaplaciansum += gammalaplacian;

            chempot +=
                2 * mv_Beta[i][j] *
                    (-12 * ci * ci * cj - 12 * ci * cj * cj + 12 * ci * cj - 4 * cj * cj * cj + 6 * cj * cj - 2 * cj) +
                gammalaplacian;
            chempot_bulk +=
                2 * mv_Beta[i][j] *
                (-12 * ci * ci * cj - 12 * ci * cj * cj + 12 * ci * cj - 4 * cj * cj * cj + 6 * cj * cj - 2 * cj);
            chempot_tension += +gammalaplacian;
        }

        const double cj = 1.0 - sumcj;

        chempot +=
            2 * mv_Beta[i][N - 1] *
                (-12 * ci * ci * cj - 12 * ci * cj * cj + 12 * ci * cj - 4 * cj * cj * cj + 6 * cj * cj - 2 * cj) +
            mv_Gamma[i][N - 1] * onlylaplaciansum;
        chempot_bulk +=
            2 * mv_Beta[i][N - 1] *
            (-12 * ci * ci * cj - 12 * ci * cj * cj + 12 * ci * cj - 4 * cj * cj * cj + 6 * cj * cj - 2 * cj);
        chempot_tension += mv_Gamma[i][N - 1] * onlylaplaciansum;

        getInstance<ChemicalPotential, N, Lattice>(i)[k] = chempot;
    }

    constexpr int i = N - 1;
    chempot = 0.0;
    chempot_bulk = 0.0;
    chempot_tension = 0.0;
    const double& ci = 1.0 - sumci;

    sumcj = 0.0;
    gammalaplaciansum = 0.0;
    onlylaplaciansum = 0.0;

    for (int j = 0; j < N - 1; j++) {
        const double& cj = getInstance<OrderParameter, N, Lattice>(j)[k];
        sumcj += cj;

        onlylaplacian = getInstance<LaplacianOrderParameter, N, Lattice>(j)[k];
        onlylaplaciansum += onlylaplacian;

        gammalaplacian = mv_Gamma[i][j] * onlylaplacian;
        gammalaplaciansum += gammalaplacian;

        if (i != j) {
            chempot +=
                2 * mv_Beta[i][j] *
                    (-12 * ci * ci * cj - 12 * ci * cj * cj + 12 * ci * cj - 4 * cj * cj * cj + 6 * cj * cj - 2 * cj) +
                gammalaplacian;
            chempot_bulk +=
                2 * mv_Beta[i][j] *
                (-12 * ci * ci * cj - 12 * ci * cj * cj + 12 * ci * cj - 4 * cj * cj * cj + 6 * cj * cj - 2 * cj);
            chempot_tension += +gammalaplacian;
        }
    }

    getInstance<ChemicalPotential, N, Lattice>(i)[k] = chempot;
}
