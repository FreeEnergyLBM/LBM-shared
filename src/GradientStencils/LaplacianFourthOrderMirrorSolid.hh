#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct LaplacianCentralFourthOrderMirrorSolid : GradientBase<One> {

    template<class TTraits, class TParameter>
    inline double compute(int direction, int k, int num = 0);

    template<class TObj>
    using GradientType = Laplacian<TObj,TObj::instances>;
    
};


template<class TTraits, class TParameter>
inline double LaplacianCentralFourthOrderMirrorSolid::compute(int direction, int k, int num){

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double laplaciansum=0;

    const static auto& param = TParameter::template get<Lattice>();

    for (int idx = 1; idx < Stencil::Q; idx++) {

            if(!this->isBoundary<Lattice>(data.getNeighbor(k, idx))&&!this->isBoundary<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx))) {

                laplaciansum +=  Stencil::Weights[idx] * 2 * (- TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx), num) + 16 * TParameter::template get<Lattice>(data.getNeighbor(k, idx), num) - 15 * TParameter::template get<Lattice>(k, num));

            }
            else if (!this->isBoundary<Lattice>(data.getNeighbor(k, idx))) {

                const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx)).NormalDirection)->second;

                double csolid = param[data.getNeighbor(data.getNeighbor(data.getNeighbor(k, idx), idx), normalq)*TParameter::instances + num];

                laplaciansum +=  Stencil::Weights[idx] * 2 * (- csolid + 16 * TParameter::template get<Lattice>(data.getNeighbor(k, idx), num) - 15 * TParameter::template get<Lattice>(k, num));

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

                laplaciansum +=  Stencil::Weights[idx] * 2 * (- csolid2 + 16 * csolid - 15 * TParameter::template get<Lattice>(k, num));

            }

            //laplaciansum +=  Stencil::Weights[idx] * 2 * (- TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx), num) + 16 * TParameter::template get<Lattice>(data.getNeighbor(k, idx), num) - 15 * TParameter::template get<Lattice>(k, num));

    }
    return 1.0 / (12.0 * Stencil::Cs2 * Lattice::DT * Lattice::DT) * laplaciansum;
}
