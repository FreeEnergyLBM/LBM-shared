#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct MixedQMirrorSolid : GradientBase<AllDirections> {

    template<class TTraits, class TParameter>
    static inline double compute(const int direction, const int k, int num = 0);

    template<class TObj>
    using GradientType = GradientMixed<TObj,TObj::instances>;
    
};

template<class TTraits, class TParameter>
inline double MixedQMirrorSolid::compute(const int direction, const int k, int num){

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (Geometry<Lattice>::getBoundaryType(k) == 1||Geometry<Lattice>::getBoundaryType(k) == 4) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    /*
    return 0.25 * (-TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q + direction], num)
                   + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num) 
                   - 3 * TParameter::template get<Lattice>(k, num) 
                   - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));
    */
    const static auto& param = TParameter::template get<Lattice>();
    /*
    if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, direction))==1)) {

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[direction]))==1)){
            const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, direction)).NormalDirection;
            const int& normalqforward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, direction)).NormalDirection)->second;
            const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;

            std::array<int8_t,TTraits::Lattice::NDIM> newdir = {};

            newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[direction]+2*(int)normal[0]*(TTraits::Stencil::Ci_x[direction]==-(int)normal[0]));
            if constexpr (TTraits::Lattice::NDIM>=2) newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[direction]+2*(int)normal[1]*(TTraits::Stencil::Ci_y[direction]==-(int)normal[1]));
            if constexpr (TTraits::Lattice::NDIM>=3) newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[direction]+2*(int)normal[2]*(TTraits::Stencil::Ci_z[direction]==-(int)normal[2]));

            const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;

            double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalqforward), num);
            double csolid2 = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalqforward]), newidx), newidx), num);
            double csolid3 = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalqbackward), num);
            //#pragma omp critical
            //{
            //std::cout<<TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num)<<" "<<TParameter::template get<Lattice>(k, num)<<" "<<csolid<<" "<<csolid2<<", "<<direction<<" "<<normalq<<" "<<newidx<<" "<<k<<std::endl;
            //}
            return 0.25 * (- (csolid2 )
                        + 5 * (csolid )
                        - 3 * TParameter::template get<Lattice>(k, num)
                        - csolid3);
        }
        else {

            const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, direction)).NormalDirection;
            const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, direction)).NormalDirection)->second;

            std::array<int8_t,TTraits::Lattice::NDIM> newdir = {};

            newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[direction]+2*(int)normal[0]*(TTraits::Stencil::Ci_x[direction]==-(int)normal[0]));
            if constexpr (TTraits::Lattice::NDIM>=2) newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[direction]+2*(int)normal[1]*(TTraits::Stencil::Ci_y[direction]==-(int)normal[1]));
            if constexpr (TTraits::Lattice::NDIM>=3) newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[direction]+2*(int)normal[2]*(TTraits::Stencil::Ci_z[direction]==-(int)normal[2]));

            const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;

            double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, direction), normalq), num);
            double csolid2 = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalq]), newidx), newidx), num);
            //#pragma omp critical
            //{
            //std::cout<<TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num)<<" "<<TParameter::template get<Lattice>(k, num)<<" "<<csolid<<" "<<csolid2<<", "<<direction<<" "<<normalq<<" "<<newidx<<" "<<k<<std::endl;
            //}
            return 0.25 * (- (csolid2 )
                        + 5 * (csolid )
                        - 3 * TParameter::template get<Lattice>(k, num)
                        - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));
        }

    }
    else if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(data.getNeighbor(k, direction), direction))==1) {
        
        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)){

            const int& normalqforward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, direction), direction)).NormalDirection)->second;
            double csolidforward = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(data.getNeighbor(k, direction), direction), normalqforward), num);
            const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;
            double csolidbackward = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalqbackward), num);

            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])==4)) {

                return 0.25 * (+ 4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                                - 3 * TParameter::template get<Lattice>(k, num)
                                - csolidbackward);
            
            }
            else return 0.25 * (- csolidforward
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - csolidbackward);
        }
        
        const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, direction), direction)).NormalDirection)->second;
        double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(data.getNeighbor(k, direction), direction), normalq), num);
        
        if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k,direction))==4) {

            return 0.25 * (+ 4 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

        }
        else return 0.25 * (- (csolid)
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    
    else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k,Stencil::Opposites[direction]))==1)) {

        const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;
        double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalq), num);

        return 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q+direction], num)
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - (csolid));

    }
    else {

        return 0.25 * (- TParameter::template get<Lattice>(data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q+direction], num)
                       + 5 * TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + direction], num)
                       - 3 * TParameter::template get<Lattice>(k, num)
                       - TParameter::template get<Lattice>(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]], num));

    }
    */

    if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])!=1)
                && (Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])!=1)) {

        return 0.25 * (- param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q+direction]*TParameter::instances + num]
                       + 5 * param[data.getNeighbors()[k * Stencil::Q + direction]*TParameter::instances + num]
                       - 3 * param[k*TParameter::instances + num]
                       - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]);

    }
    else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, direction))==1)) {

        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[direction]))==1)){
            const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, direction)).NormalDirection;
            const int& normalqforward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, direction)).NormalDirection)->second;
            const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;

            std::array<int8_t,TTraits::Lattice::NDIM> newdir = {};

            newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[direction]+2*(int)normal[0]*(TTraits::Stencil::Ci_x[direction]==-(int)normal[0]));
            if constexpr (TTraits::Lattice::NDIM>=2) newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[direction]+2*(int)normal[1]*(TTraits::Stencil::Ci_y[direction]==-(int)normal[1]));
            if constexpr (TTraits::Lattice::NDIM>=3) newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[direction]+2*(int)normal[2]*(TTraits::Stencil::Ci_z[direction]==-(int)normal[2]));

            const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;

            double csolid = param[data.getNeighbor(data.getNeighbor(k, direction), normalqforward)*TParameter::instances + num];
            double csolid2 = param[data.getNeighbor(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalqforward]), newidx), newidx)*TParameter::instances + num];
            double csolid3 = param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalqbackward)*TParameter::instances + num];
            //#pragma omp critical
            //{
            //std::cout<<param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]<<" "<<param[k*TParameter::instances + num]<<" "<<csolid<<" "<<csolid2<<", "<<direction<<" "<<normalq<<" "<<newidx<<" "<<k<<std::endl;
            //}
            return 0.25 * (- (csolid2 )
                        + 5 * (csolid )
                        - 3 * param[k*TParameter::instances + num]
                        - csolid3);
        }
        else {

            const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, direction)).NormalDirection;
            const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, direction)).NormalDirection)->second;

            std::array<int8_t,TTraits::Lattice::NDIM> newdir = {};

            newdir[0] = (int8_t)(TTraits::Stencil::Ci_x[direction]+2*(int)normal[0]*(TTraits::Stencil::Ci_x[direction]==-(int)normal[0]));
            if constexpr (TTraits::Lattice::NDIM>=2) newdir[1] = (int8_t)(TTraits::Stencil::Ci_y[direction]+2*(int)normal[1]*(TTraits::Stencil::Ci_y[direction]==-(int)normal[1]));
            if constexpr (TTraits::Lattice::NDIM>=3) newdir[2] = (int8_t)(TTraits::Stencil::Ci_z[direction]+2*(int)normal[2]*(TTraits::Stencil::Ci_z[direction]==-(int)normal[2]));

            const int& newidx = TTraits::Stencil::QMap.find(newdir)->second;

            double csolid = param[data.getNeighbor(data.getNeighbor(k, direction), normalq)*TParameter::instances + num];
            double csolid2 = param[data.getNeighbor(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[normalq]), newidx), newidx)*TParameter::instances + num];
            //#pragma omp critical
            //{
            //std::cout<<param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]<<" "<<param[k*TParameter::instances + num]<<" "<<csolid<<" "<<csolid2<<", "<<direction<<" "<<normalq<<" "<<newidx<<" "<<k<<std::endl;
            //}
            return 0.25 * (- (csolid2 )
                        + 5 * (csolid )
                        - 3 * param[k*TParameter::instances + num]
                        - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]);
        }

    }
    else if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(data.getNeighbor(k, direction), direction))==1) {
        
        if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]])==1)){

            const int& normalqforward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, direction), direction)).NormalDirection)->second;
            double csolidforward = param[data.getNeighbor(data.getNeighbor(data.getNeighbor(k, direction), direction), normalqforward)*TParameter::instances + num];
            const int& normalqbackward = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;
            double csolidbackward = param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalqbackward)*TParameter::instances + num];

            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbors()[k * Stencil::Q + direction])==4)) {

                return 0.25 * (+ 4 * param[data.getNeighbors()[k * Stencil::Q + direction]*TParameter::instances + num]
                                - 3 * param[k*TParameter::instances + num]
                                - csolidbackward);
            
            }
            else return 0.25 * (- csolidforward
                       + 5 * param[data.getNeighbors()[k * Stencil::Q + direction]*TParameter::instances + num]
                       - 3 * param[k*TParameter::instances + num]
                       - csolidbackward);
        }
        
        const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, direction), direction)).NormalDirection)->second;
        double csolid = param[data.getNeighbor(data.getNeighbor(data.getNeighbor(k, direction), direction), normalq)*TParameter::instances + num];
        
        if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k,direction))==4) {

            return 0.25 * (+ 4 * param[data.getNeighbors()[k * Stencil::Q + direction]*TParameter::instances + num]
                       - 3 * param[k*TParameter::instances + num]
                       - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]);

        }
        else return 0.25 * (- (csolid)
                       + 5 * param[data.getNeighbors()[k * Stencil::Q + direction]*TParameter::instances + num]
                       - 3 * param[k*TParameter::instances + num]
                       - param[data.getNeighbors()[k * Stencil::Q + Stencil::Opposites[direction]]*TParameter::instances + num]);

    }
    
    else {

        const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, Stencil::Opposites[direction])).NormalDirection)->second;
        double csolid = param[data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[direction]), normalq)*TParameter::instances + num];

        return 0.25 * (- param[data.getNeighbors()[data.getNeighbors()[k * Stencil::Q + direction] * Stencil::Q+direction]*TParameter::instances + num]
                       + 5 * param[data.getNeighbors()[k * Stencil::Q + direction]*TParameter::instances + num]
                       - 3 * param[k*TParameter::instances + num]
                       - (csolid));

    }
    
}