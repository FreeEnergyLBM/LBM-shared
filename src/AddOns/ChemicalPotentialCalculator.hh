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
        Lattice::communicate(ChemicalPotential<3>::getInstance<Lattice>());
    }

    std::array<double, 3> mSigma = {};
    std::array<double, 3> mGamma = {};
    double mGammaT = 0;

    double mOmega = 0.0001;
    double mLambda = 0.00;
    double mInterfaceWidth = 4;
};

template <class TTraits>
inline void ChemicalPotentialCalculatorTernaryLee::compute(int k) {
    using Lattice = typename TTraits::Lattice;

    if (Geometry<Lattice>::isBulkSolid(k)) return;

    const double C1 = OrderParameter<2>::get<Lattice>(k, 0);  // C1, C2, C3
    const double C2 = OrderParameter<2>::get<Lattice>(k, 1);
    const double C3 = 1 - C1 - C2;

    double dEd12 = mGamma[0] * C1 * (1 - 3 * C1 + 2 * C1 * C1) - mGamma[1] * C2 * (1 - 3 * C2 + 2 * C2 * C2) +
                   2 * mLambda * C3 * C3 * (C2 * C2 * C1 - C1 * C1 * C2) + (C1 < 0) * 2 * mOmega * C1 -
                   (C2 < 0) * 2 * mOmega * C2;

    double dEd13 = mGamma[0] * C1 * (1 - 3 * C1 + 2 * C1 * C1) - mGamma[2] * C3 * (1 - 3 * C3 + 2 * C3 * C3) +
                   2 * mLambda * C2 * C2 * (C3 * C3 * C1 - C1 * C1 * C3) + (C1 < 0) * 2 * mOmega * C1 -
                   (C3 < 0) * 2 * mOmega * C3;

    double dEd23 = mGamma[1] * C2 * (1 - 3 * C2 + 2 * C2 * C2) - mGamma[2] * C3 * (1 - 3 * C3 + 2 * C3 * C3) +
                   2 * mLambda * C1 * C1 * (C3 * C3 * C2 - C2 * C2 * C3) + (C2 < 0) * 2 * mOmega * C2 -
                   (C3 < 0) * 2 * mOmega * C3;

    ChemicalPotential<3>::get<Lattice>(k, 0) =
        (4.0 * mGammaT / mInterfaceWidth * (1 / mGamma[1] * dEd12 + 1 / mGamma[2] * dEd13) -
         3.0 / 4.0 * mInterfaceWidth * mGamma[0] * (LaplacianOrderParameter<2>::get<Lattice>(k, 0)));

    ChemicalPotential<3>::get<Lattice>(k, 1) =
        (4.0 * mGammaT / mInterfaceWidth * (-1 / mGamma[0] * dEd12 + 1 / mGamma[2] * dEd23) -
         3.0 / 4.0 * mInterfaceWidth * mGamma[1] * (LaplacianOrderParameter<2>::get<Lattice>(k, 1)));

    ChemicalPotential<3>::get<Lattice>(k, 2) =
        (4.0 * mGammaT / mInterfaceWidth * (-1 / mGamma[0] * dEd13 - 1 / mGamma[1] * dEd23) -
         3.0 / 4.0 * mInterfaceWidth * mGamma[2] *
             (-LaplacianOrderParameter<2>::get<Lattice>(k, 0) - LaplacianOrderParameter<2>::get<Lattice>(k, 1)));
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

class ChemicalPotentialCalculatorNComponent : public AddOnBase {
   public:
    ChemicalPotentialCalculatorNComponent() : mv_Beta({{0}}), mv_Gamma({{0}}) {}

    double** ma_Gamma;

    double** ma_Beta;

    std::vector<std::vector<double>> mv_Beta;
    std::vector<std::vector<double>> mv_Gamma;

    double mD = 5;

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
};

template <class TTraits>
inline void ChemicalPotentialCalculatorNComponent::compute(const int k) {
    if (Geometry<typename TTraits::Lattice>::isBulkSolid(k)) return;

    [[maybe_unused]] static bool isvalid = checkValid<TTraits::NumberOfComponents>();
    mIsValid = isvalid;

    double chempot = 0.0;
    double chempot_bulk = 0.0;
    double chempot_tension = 0.0;
    double sumci = 0.0;
    double sumcj = 0.0;

    double onlylaplacian = 0.0;
    double gammalaplacian = 0.0;
    double onlylaplaciansum = 0.0;
    double gammalaplaciansum = 0.0;

    for (int i = 0; i < TTraits::NumberOfComponents - 1; i++) {
        chempot = 0.0;
        chempot_bulk = 0.0;
        chempot_tension = 0.0;
        const double& ci = OrderParameter<TTraits::NumberOfComponents - 1>::template get<class TTraits::Lattice>(k, i);
        sumci += ci;

        sumcj = 0.0;
        gammalaplaciansum = 0.0;
        onlylaplaciansum = 0.0;

        for (int j = 0; j < TTraits::NumberOfComponents - 1; j++) {
            const double& cj =
                OrderParameter<TTraits::NumberOfComponents - 1>::template get<class TTraits::Lattice>(k, j);
            sumcj += cj;

            onlylaplacian =
                LaplacianOrderParameter<TTraits::NumberOfComponents - 1>::template get<class TTraits::Lattice>(k, j);

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
            2 * mv_Beta[i][TTraits::NumberOfComponents - 1] *
                (-12 * ci * ci * cj - 12 * ci * cj * cj + 12 * ci * cj - 4 * cj * cj * cj + 6 * cj * cj - 2 * cj) +
            mv_Gamma[i][TTraits::NumberOfComponents - 1] * onlylaplaciansum;
        chempot_bulk +=
            2 * mv_Beta[i][TTraits::NumberOfComponents - 1] *
            (-12 * ci * ci * cj - 12 * ci * cj * cj + 12 * ci * cj - 4 * cj * cj * cj + 6 * cj * cj - 2 * cj);
        chempot_tension += mv_Gamma[i][TTraits::NumberOfComponents - 1] * onlylaplaciansum;

        ChemicalPotential<TTraits::NumberOfComponents>::template get<class TTraits::Lattice>(k, i) = chempot;
    }

    int i = TTraits::NumberOfComponents - 1;
    chempot = 0.0;
    chempot_bulk = 0.0;
    chempot_tension = 0.0;
    const double& ci = 1.0 - sumci;

    sumcj = 0.0;
    gammalaplaciansum = 0.0;
    onlylaplaciansum = 0.0;

    for (int j = 0; j < TTraits::NumberOfComponents - 1; j++) {
        const double& cj = OrderParameter<TTraits::NumberOfComponents - 1>::template get<class TTraits::Lattice>(k, j);
        sumcj += cj;

        onlylaplacian =
            LaplacianOrderParameter<TTraits::NumberOfComponents - 1>::template get<class TTraits::Lattice>(k, j);
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

    ChemicalPotential<TTraits::NumberOfComponents>::template get<class TTraits::Lattice>(k, i) = chempot;
}
