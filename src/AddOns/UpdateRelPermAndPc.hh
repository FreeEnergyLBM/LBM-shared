#pragma once
#include "../Geometry.hh"
#include "../Lattice.hh"
#include "../Parameters.hh"
#include "../Service.hh"
#include "AddOnBase.hh"

/** This class uses the Brooks-Corey model to update the relative permeability and capillary pressure
 * based on the saturation field. The capillary pressure and relative permeability equations can be found
 * in 10.1029/2005WR004482.
 *
 * For cases where the porosity and permeability fields are variable, the Leverett J-function is used to
 * scale the capillary pressure.
 */

class UpdateRelPermAndPc : public AddOnBase {
   public:
    // Functions to set the irreducible wetting and non-wetting saturations
    inline void setSwIr(double SwIr) { Sw_ir = SwIr; }
    inline void setSnwIr(double SnwIr) { Snw_ir = SnwIr; }

    // Functions to set the relative permeability function exponents
    inline void setN(double N) { n = N; }
    inline void setM(double M) { m = M; }

    // Functions to set the capillary pressure parameters
    inline void setP0(double P0) { P_0 = P0; }
    inline void setPcInf(double PcInf) { P_c_inf = PcInf; }
    inline void setLambda(double Lambda) { lambda = Lambda; }

    // Functions to set the Leverett J-function parameters
    inline void setK0(double K0) { K_0 = K0; }
    inline void setPhi0(double Phi0) { phi_0 = Phi0; }

    template <class TTraits>
    inline void compute(int k);

   private:
    inline double computeCapillaryPressure(double SwStar, double LJFunctionPrefactor);
    inline double computeRelativePermeability(double SwStar);

    // Irreducible wetting and non-wetting saturations
    double Sw_ir = 0.0;
    double Snw_ir = 0.0;

    // Leverett J-function parameters
    double K_0 = 1.0;
    double phi_0 = 1.0;

    // Capillary pressure parameters
    double P_0 = 1.0E-2;
    double P_c_inf = 10 * P_0;
    double lambda = 2.0;

    // Default values for the relative permeability exponents (Brooks-Corey model)
    double n = (2 + 3 * lambda) / lambda;
    double m = (2 + 3 * lambda) / lambda;
};

template <class TTraits>
inline void UpdateRelPermAndPc::compute(int k) {
    using Lattice = typename TTraits::Lattice;

    if (Geometry<Lattice>::isBulkSolid(k)) return;

    double SwStar = (Saturation<>::get<Lattice>(k) - Sw_ir) / (1 - Sw_ir - Snw_ir);
    if (SwStar < 0) SwStar = 0;
    if (SwStar > 1) SwStar = 1;

    // Update Capillary Pressure
    double porosityRatio = Porosity<>::get<Lattice>(k) / phi_0;
    double permeabilityRatio = Permeability<>::get<Lattice>(k) / K_0;

    double LJFunctionPrefactor = sqrt(porosityRatio / permeabilityRatio);
    double Pc = computeCapillaryPressure(SwStar, LJFunctionPrefactor);
    CapillaryPressure<>::get<Lattice>(k) = Pc;

    // Update Relative Permeability
    double relPerm = computeRelativePermeability(SwStar);
    WettingRelativePermeability<>::get<Lattice>(k) = relPerm;
}

inline double UpdateRelPermAndPc::computeCapillaryPressure(double SwStar, double LJFunctionPrefactor) {
    // If condition is added to avoid division by zero and very high capillary pressures
    if (SwStar < 0.01) return P_c_inf;

    // Scale the capillary pressure using the Leverett J-function
    return LJFunctionPrefactor * P_0 * pow(SwStar, -1 / lambda);
}

inline double UpdateRelPermAndPc::computeRelativePermeability(double SwStar) { return pow(SwStar, n); }
