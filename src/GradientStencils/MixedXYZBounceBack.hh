#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedXYZBounceBack : GradientBase<GradientMixed,Cartesian> {

    template<class TTraits, class TParameter>
    inline double compute(const int direction, const int k, int num = 0);
    
};

template<class TTraits, class TParameter>
inline double MixedXYZBounceBack::compute(const int direction, const int k, int num){

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    if (this->isBoundary<Lattice>(k)) return 0;

    static DataType& data = DataType::getInstance();

    double gradientsum=0;
    const static auto& param = TParameter::template get<Lattice>();
    for (int idx = 1; idx < Stencil::Q; idx++) {

        if ((!this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + idx]))&& (!this->isBoundary<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx]))
                && (!this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {
                
            
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx]
                                                                                                                                                * Stencil::Q+  idx]*TParameter::instances + num]
                                                                                       + 5 * param[data.getNeighbor(k, idx)*TParameter::instances + num]
                                                                                       - 3 * param[k*TParameter::instances + num]
                                                                                       - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]*TParameter::instances + num]);
            /*
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx]
                                                                                                                                                * Stencil::Q+  idx], num)
                                                                                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                                                                       - 3 * TParameter::template get<Lattice>(k, num)
                                                                                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));
            */

        }
        else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + idx]))) {
            /*
            if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]])!=1)) gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (2 *  TParameter::template get<Lattice>(k, num)
                                                                                       - 2 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));
            */
            if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]])!=1)) gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (2 *  param[k *TParameter::instances + num]
                                                                                       - 2 * param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]] *TParameter::instances + num]);

        }
        else if ((this->isBoundary<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx]))) {
            /*
            if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]])==1)) {
                
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num) 
                                                                                       - 4 * TParameter::template get<Lattice>(k, num));

            }
            else gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num) 
                                                                                       - 3 * TParameter::template get<Lattice>(k, num)
                                                                                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));
            */
            if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {
                
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (4 * param[data.getNeighbors()[k * Stencil::Q + idx] *TParameter::instances + num] 
                                                                                       - 4 * param[k *TParameter::instances + num]);

            }
            else gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (4 * param[data.getNeighbors()[k * Stencil::Q + idx] *TParameter::instances + num] 
                                                                                       - 3 * param[k *TParameter::instances + num]
                                                                                       - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]] *TParameter::instances + num]);
        }
        else {
            /*
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx]
                                                                                                                                                * Stencil::Q + idx], num)
                                                                                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                                                                       - 4 * TParameter::template get<Lattice>(k, num));
            */
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx]
                                                                                                                                                * Stencil::Q + idx]*TParameter::instances + num]
                                                                                       + 5 * param[data.getNeighbors()[k * Stencil::Q + idx]*TParameter::instances + num]
                                                                                       - 4 * param[k*TParameter::instances + num]);
        }

    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
        
}