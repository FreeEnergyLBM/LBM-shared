#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct LaplacianCentralMirrorSolid : GradientBase<Laplacian,One> { //FIX

    template<class TTraits, class TParameter>
    inline double compute(int direction, int k, int num = 0);
    
};

template<class TTraits, class TParameter>
inline double LaplacianCentralMirrorSolid::compute(const int direction, const int k, int num){
   
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    //if (Geometry<Lattice>::getBoundaryType(k) == 4) return 0;
    if (this->isBoundary<Lattice>(k)) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double laplaciansum=0;
    const static auto& param = TParameter::template get<Lattice>();
    /*
    for (int idx = 1; idx <Stencil::Q; idx++) {

        if(Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))!=1) {

            laplaciansum +=  Stencil::Weights[idx] * 2 * (TParameter::template get<Lattice>(data.getNeighbor(k, idx), num) - TParameter::template get<Lattice>(k, num));

        }
        else {

            const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;

            double csolid = TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), normalq), num);

            laplaciansum +=  Stencil::Weights[idx] * 2 * (csolid - TParameter::template get<Lattice>(k, num));

        }

    }
    */
    for (int idx = 1; idx <Stencil::Q; idx++) {

        //if (Stencil::Ci_xyz(1)[idx]!=0) continue;

        if(Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))!=1) {

            laplaciansum +=  Stencil::Weights[idx] * 2 * (param[data.getNeighbor(k, idx)*TParameter::instances + num] - param[k*TParameter::instances + num]);

        }
        else {

            const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(data.getNeighbor(k, idx)).NormalDirection)->second;

            double csolid = param[data.getNeighbor(data.getNeighbor(k, idx), normalq)*TParameter::instances + num];

            laplaciansum +=  Stencil::Weights[idx] * 2 * (csolid - param[k*TParameter::instances + num]);

        }

    }
    return 1.0 / (Stencil::Cs2 * Lattice::DT * Lattice::DT) * laplaciansum;

}
