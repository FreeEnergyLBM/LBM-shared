#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct LaplacianCentralWetting : WettingGradient<One> {

    template<class TTraits, class TParameter>
    inline double compute(int direction, int k, int num = 0);

    template<class TObj>
    using GradientType = Laplacian<TObj,TObj::instances>;
    
};


template<class TTraits, class TParameter>
inline double LaplacianCentralWetting::compute(int direction, int k, int num){

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (Geometry<Lattice>::getBoundaryType(k) == 4) return 0;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double laplaciansum=0;

    for (int idx = 1; idx <Stencil::Q; idx++) {

        if(Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))!=1) {

            laplaciansum +=  Stencil::Weights[idx] * 2 * (TParameter::template get<Lattice>(data.getNeighbor(k, idx), num) - TParameter::template get<Lattice>(k, num));

        }
        else {

            double csolid = TParameter::template get<Lattice>(k, num);//TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), normalq), num);
            //#pragma omp critical
            //{
            //if(Geometry<Lattice>::getBoundaryType(k)!=1) std::cout<<normalq<<" "<<k<<" "<<idx<<" "<<Geometry<Lattice>::getBoundaryType(k)<<std::endl;
            //}
            laplaciansum +=  Stencil::Weights[idx] * 2 * ((csolid - 0.5 * this->mPrefactor * (csolid - pow(csolid, 2))) - TParameter::template get<Lattice>(k, num));

        }

    }
    return 1.0 / (Stencil::Cs2 * Lattice::DT * Lattice::DT) * laplaciansum;
}
