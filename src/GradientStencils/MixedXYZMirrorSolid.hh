#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedXYZMirrorSolid : GradientBase<Cartesian> {

    template<class TTraits, class TParameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class TObj>
    using GradientType = GradientMixed<TObj,TObj::instances>;
    
};

template<class TTraits, class TParameter>
inline double MixedXYZMirrorSolid::compute(const int direction, const int k, int num){

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (Geometry<Lattice>::getBoundaryType(k) == 4) return 0;
    
    using DataType = Data_Base<Lattice, Stencil>;

    static DataType& data = DataType::getInstance();
    const static auto& param = TParameter::template get<Lattice>();
    double gradientsum=0;
    /*
    for (int idx = 1; idx < Stencil::Q; idx++) {
        //if (idx == 2) std::cerr<<" s "<<data.getNeighbors().size()<<" n1 "<<data.getNeighbors()[k * Stencil::Q + idx]<<" n2 "<<data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx]<<std::endl;
        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + idx])==1)) {
            
            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]])==1)) {
                
                const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection;
                const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;
                const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[idx])).NormalDirection)->second;

                std::array<int8_t,TTraits::Lattice::NDIM> newdir = {};

                newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[idx]+2*(int)normal[0]*(TTraits::Stencil::Ci_x[idx]==-(int)normal[0]));
                if constexpr (TTraits::Lattice::NDIM>=2) newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[idx]+2*(int)normal[1]*(TTraits::Stencil::Ci_y[idx]==-(int)normal[1]));
                if constexpr (TTraits::Lattice::NDIM>=3) newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[idx]+2*(int)normal[2]*(TTraits::Stencil::Ci_z[idx]==-(int)normal[2]));

                const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;
                
                double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), normalq), num);
                double csolid2 = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalq]), newidx), newidx), num);
                double csolid3 = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), normalqbackward), num);
                //const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;
                //double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), normalq), num);
                
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- (csolid2)
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
                
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- (csolid2)
                                                                                            + 5 * (csolid)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));

            }
            
        }
        else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx])==1)) {
            
            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]])==1)) {
                
                const int& normalqforward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx)).NormalDirection)->second;
                double csolidforward = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(data.getNeighbor(k, idx), idx), normalqforward), num);
                const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[idx])).NormalDirection)->second;
                double csolidbackward = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), normalqbackward), num);

                if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + idx])==4)) {

                    return gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (+ 4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                - 3 * TParameter::template get<Lattice>(k, num)
                                - (csolidbackward));
            
                }
                else gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- (csolidforward)
                                                                                            + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - csolidbackward);
                
            }
            else {

                const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx)).NormalDirection)->second;
                double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(data.getNeighbor(k, idx), idx), normalq), num);
                
                if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + idx])==4)) {

                    gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));

                }
                else gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- (csolid)
                                                                                            + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                                                                            - 3 * TParameter::template get<Lattice>(k, num)
                                                                                            - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));
                                                                                    
            }
            
        }
        else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]])==1)) {
            
            const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[idx])).NormalDirection)->second;
            double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), normalq), num);

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q+idx], num)
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - (csolid));

        }
        else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + idx])!=1)
                || (Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]])!=1)) {
            
            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx]
                                                                                                                                                * Stencil::Q+  idx], num)
                                                                                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + idx], num)
                                                                                       - 3 * TParameter::template get<Lattice>(k, num)
                                                                                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]], num));

        }

    }
    */

    for (int idx = 1; idx < Stencil::Q; idx++) {

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + idx])!=1)
                && (Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]])!=1)) {

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx]
                                                                                                                                                * Stencil::Q+  idx]*TParameter::instances + num]
                                                                                       + 5 * param[data.getNeighbors()[k * Stencil::Q + idx]*TParameter::instances + num]
                                                                                       - 3 * param[k*TParameter::instances + num]
                                                                                       - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]*TParameter::instances + num]);

        }
        else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + idx])==1)) {

            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]])==1)) {

                const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection;
                const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;
                const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[idx])).NormalDirection)->second;

                std::array<int8_t,TTraits::Lattice::NDIM> newdir = {};

                newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[idx]+2*(int)normal[0]*(TTraits::Stencil::Ci_x[idx]==-(int)normal[0]));
                if constexpr (TTraits::Lattice::NDIM>=2) newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[idx]+2*(int)normal[1]*(TTraits::Stencil::Ci_y[idx]==-(int)normal[1]));
                if constexpr (TTraits::Lattice::NDIM>=3) newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[idx]+2*(int)normal[2]*(TTraits::Stencil::Ci_z[idx]==-(int)normal[2]));

                const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;

                double csolid = param[data.getNeighbor(data.getNeighbor(k, idx), normalq)*TParameter::instances + num];
                double csolid2 = param[data.getNeighbor(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalq]), newidx), newidx)*TParameter::instances + num];
                double csolid3 = param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), normalqbackward)*TParameter::instances + num];
                //const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;
                //double csolid = param[data.getNeighbor(data.getNeighbor(k, idx), normalq)*TParameter::instances + num];
                
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- (csolid2)
                                                                                            + 5 * (csolid)
                                                                                            - 3 * param[k*TParameter::instances + num]
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

                double csolid = param[data.getNeighbor(data.getNeighbor(k, idx), normalq)*TParameter::instances + num];
                double csolid2 = param[data.getNeighbor(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalq]), newidx), newidx)*TParameter::instances + num];

                //const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;
                //double csolid = param[data.getNeighbor(data.getNeighbor(k, idx), normalq)*TParameter::instances + num];
                
                gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- (csolid2)
                                                                                            + 5 * (csolid)
                                                                                            - 3 * param[k*TParameter::instances + num]
                                                                                            - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]*TParameter::instances + num]);

            }

        }
        else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q + idx])==1)) {

            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]])==1)) {

                const int& normalqforward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx)).NormalDirection)->second;
                double csolidforward = param[data.getNeighbor(data.getNeighbor(data.getNeighbor(k, idx), idx), normalqforward)*TParameter::instances + num];
                const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[idx])).NormalDirection)->second;
                double csolidbackward = param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), normalqbackward)*TParameter::instances + num];

                if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + idx])==4)) {

                    return gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (+ 4 * param[data.getNeighbors()[k * Stencil::Q + idx]*TParameter::instances + num]
                                - 3 * param[k*TParameter::instances + num]
                                - (csolidbackward));
            
                }
                else gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- (csolidforward)
                                                                                            + 5 * param[data.getNeighbors()[k * Stencil::Q + idx]*TParameter::instances + num]
                                                                                            - 3 * param[k*TParameter::instances + num]
                                                                                            - csolidbackward);

            }
            else {

                const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx)).NormalDirection)->second;
                double csolid = param[data.getNeighbor(data.getNeighbor(data.getNeighbor(k, idx), idx), normalq)*TParameter::instances + num];
                
                if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + idx])==4)) {

                    gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (4 * param[data.getNeighbors()[k * Stencil::Q + idx]*TParameter::instances + num]
                                                                                            - 3 * param[k*TParameter::instances + num]
                                                                                            - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]*TParameter::instances + num]);

                }
                else gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- (csolid)
                                                                                            + 5 * param[data.getNeighbors()[k * Stencil::Q + idx]*TParameter::instances + num]
                                                                                            - 3 * param[k*TParameter::instances + num]
                                                                                            - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[idx]]*TParameter::instances + num]);
                                                                                    
            }

        }
        else {

            const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[idx])).NormalDirection)->second;
            double csolid = param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), normalq)*TParameter::instances + num];

            gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * 0.25 * (- param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + idx] * Stencil::Q+idx]*TParameter::instances + num]
                       + 5 * param[data.getNeighbors()[k * Stencil::Q + idx]*TParameter::instances + num]
                       - 3 * param[k*TParameter::instances + num]
                       - (csolid));

        }

    }
    
    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;
    
}