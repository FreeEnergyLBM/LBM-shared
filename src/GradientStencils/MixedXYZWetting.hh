#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedXYZWetting : WettingGradient<Cartesian> {

    template<class TTraits, class TParameter>
    inline double compute(const int direction, const int k, int num = 0);

    template<class TObj>
    using GradientType = GradientMixed<TObj,TObj::instances>;

};

template<class TTraits, class TParameter>
inline double MixedXYZWetting::compute(const int direction, const int k, int num){
    
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (Geometry<Lattice>::getBoundaryType(k) == 4) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum=0;

    for (int idx = 1; idx < Stencil::Q; idx++) {

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])==1)) {

            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])!=1)) {

                double csolid = TParameter::template get<Lattice>(k, num);
                
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- (csolid - 1.5 * this->mPrefactor * (csolid - pow(csolid, 2)))
                                                                                            + 5 * (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2)))
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

            }

        }
        else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction])==1)) {

            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)) {

                double csolid = TParameter::template get<Lattice>(k, num);

                if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])==4)) {

                    return gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * (+ 4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                                - 3 * TParameter::template get<Lattice>(k, num)
                                - (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))));
            
                }
                else gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2)))
                                                                                            + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))));

            }
            else {

                double csolid = TParameter::template get<Lattice>(k, num);
                
                if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])==4)) {

                    gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

                }
                else gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2)))
                                                                                            + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));
                                                                                    
            }

        }
        else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)) {

            double csolid = TParameter::template get<Lattice>(k, num);

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q+direction], num)
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))));

        }
        else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])!=1)
                || (Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])!=1)) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction]
                                                                                                                                                * Stencil::Q+  direction], num)
                                                                                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                                                                                       - 3 * TParameter::template get<Lattice>(k, num)
                                                                                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

        }

    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;

}