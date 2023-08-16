#pragma once
#include "../Parameters.hh"
#include "AddOnBase.hh"
#include<iostream>
#include<math.h>


class ConstantPressureGradientBoundary : public AddOnBase {
    public:

        template<class TTraits>
        inline void compute(int k);

        int mBoundaryLabel=4;

};

template<class TTraits>
inline void ConstantPressureGradientBoundary::compute(int k) {

    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (Geometry<Lattice>::getBoundaryType(k)!=mBoundaryLabel) return;

    using data = Data_Base<Lattice, Stencil>;

    std::vector<int>& neighbors = data::getInstance().getNeighbors();
    std::vector<int8_t>& normal = BoundaryLabels<>::get<typename TTraits::Lattice>(k).NormalDirection;
    int idx = Stencil::QMap.find(normal)->second;
    //std::cout<<idx<<" "<<(int)(normal)[0]<<" "<<(int)(normal)[1]<<std::endl;
    Pressure<>::get<Lattice>(k) = (4*Pressure<>::get<Lattice>(neighbors[k * Stencil::Q+idx]) - Pressure<>::get<Lattice>(neighbors[neighbors[k * Stencil::Q+idx] * Stencil::Q+idx])) / 3.0;
    
}