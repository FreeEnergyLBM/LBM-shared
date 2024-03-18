#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedXYZWetting : WettingGradient<Cartesian> {

    template<class TTraits, class TParameter>
    inline double compute(const int direction, const int k, int num = 0);

    template<class TObj>
    using GradientType = GradientMixed<TObj,TObj::instances>;

};
/*
template<class TTraits, class TParameter>
inline double MixedXYZWetting::compute(const int direction, const int k, int num){
    
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    //if (this->isBoundary<Lattice>(k) == 4) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum=0;

    for (int idx = 1; idx < Stencil::Q; idx++) {

        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction])==1)) {

            if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])!=1)) {

                double csolid = TParameter::template get<Lattice>(k, num);
                
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- (csolid - 1.5 * this->mPrefactor * (csolid - pow(csolid, 2)))
                                                                                            + 5 * (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2)))
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

            }

        }
        else if ((this->isBoundary<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction])==1)) {

            if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)) {

                double csolid = TParameter::template get<Lattice>(k, num);

                if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction])==4)) {

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
                
                if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction])==4)) {

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
        else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)) {

            double csolid = TParameter::template get<Lattice>(k, num);

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q+direction], num)
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - (csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))));

        }
        else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + direction])!=1)
                || (this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])!=1)) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction]
                                                                                                                                                * Stencil::Q+  direction], num)
                                                                                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                                                                                       - 3 * TParameter::template get<Lattice>(k, num)
                                                                                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

        }

    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;

}
*/


template<class TTraits, class TParameter>
inline double MixedXYZWetting::compute(const int direction, const int k, int num){
        
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (this->isBoundary<Lattice>(k)) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum=0;

    for (int idx = 1; idx < Stencil::Q; idx++) {

        if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + idx]))) {

            if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {

                const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection;
                const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;
                const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;

                std::array<int8_t,TTraits::Lattice::NDIM> newdir = {};

                newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[idx]+2*(int)normal[0]*(TTraits::Stencil::Ci_x[idx]==-(int)normal[0]));
                if constexpr (TTraits::Lattice::NDIM>=2) newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[idx]+2*(int)normal[1]*(TTraits::Stencil::Ci_y[idx]==-(int)normal[1]));
                if constexpr (TTraits::Lattice::NDIM>=3) newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[idx]+2*(int)normal[2]*(TTraits::Stencil::Ci_z[idx]==-(int)normal[2]));

                const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;

                double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), normalq), num);
                double csolid2 = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalq]), newidx), newidx), num);
                double csolid3 = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalqbackward), num);
                //const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;
                //double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), normalq), num);
                
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- (csolid2)
                                                                                            + 5 * (csolid)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - csolid3);

            }
            else {

                const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection;
                const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;

                std::array<int8_t,TTraits::Lattice::NDIM> newdir = {};

                newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[idx]+2*(int)normal[0]*(TTraits::Stencil::Ci_x[idx]==-(int)normal[0]));
                if constexpr (TTraits::Lattice::NDIM>=2) newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[idx]+2*(int)normal[1]*(TTraits::Stencil::Ci_y[idx]==-(int)normal[1]));
                if constexpr (TTraits::Lattice::NDIM>=3) newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[idx]+2*(int)normal[2]*(TTraits::Stencil::Ci_z[idx]==-(int)normal[2]));

                const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;

                double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), normalq), num);
                double csolid2 = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalq]), newidx), newidx), num);

                //const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;
                //double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), normalq), num);
                
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- (csolid2)
                                                                                            + 5 * (csolid)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));

            }

        }
        else if ((this->isBoundary<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx]))) {

            if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {

                const int& normalqforward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx)).NormalDirection)->second;
                double csolidforward = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(data.getNeighbor(k, idx), idx), normalqforward), num);
                const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[idx])).NormalDirection)->second;
                double csolidbackward = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), normalqbackward), num);

                /*if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + idx])==4)) {

                    return gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * (+ 4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                - 3 * TParameter::template get<Lattice>(k, num)
                                - (csolidbackward));
            
                }*/
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- (csolidforward)
                                                                                            + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - csolidbackward);

            }
            else {

                const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx)).NormalDirection)->second;
                double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(data.getNeighbor(k, idx), idx), normalq), num);
                
                /*if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + idx])==4)) {

                    gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));

                }*/
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- (csolid)
                                                                                            + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));
                                                                                    
            }

        }
        else if ((this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {

            const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[idx])).NormalDirection)->second;
            double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), normalq), num);

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q+idx], num)
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - (csolid));

        }
        else if ((!this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + idx]))
                || (!this->isBoundary<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]))) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(idx)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx]
                                                                                                                                                * Stencil::Q+  idx], num)
                                                                                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                                                                       - 3 * TParameter::template get<Lattice>(k, num)
                                                                                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));

        }

    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
        
}