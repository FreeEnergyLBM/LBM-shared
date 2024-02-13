#pragma once
#include "BoundaryBase.hh"
#include "../Geometry.hh"
#include <iostream>

template<class TParam0thMoment>
class Neumann : public BoundaryBase {
    public:

        Neumann() { this->setNodeID(4, true); } // TMP: Default NodeID warning

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        inline void setNormalGradient(int normalgradient) {mNormalGradient=normalgradient;}

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);

    private:

        double mNormalGradient = 0;

};

template<class TParam0thMoment>
template<class TTraits, class TDistributionType>
inline void Neumann<TParam0thMoment>::compute(TDistributionType& distribution, int k) {

    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;

    const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;
    const int& normalq = TTraits::Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection)->second;
    
    const double oldval = TParam0thMoment::template get<TTraits::Lattice>(k);

    double distsum = 0;
    double weightsum = 0;
    double prefactorsum = 0;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {
        double cidotnormal = TTraits::Stencil::Ci_x[idx] * normal[0] + TTraits::Stencil::Ci_y[idx] * normal[1] * (TTraits::Lattice::NDIM>1) + TTraits::Stencil::Ci_z[idx] * normal[2] * (TTraits::Lattice::NDIM>2);
        
        if(cidotnormal > 0 ) {
            
            weightsum += TTraits::Stencil::Weights[idx];

        }
        else if ( idx>0 && Geometry<typename TTraits::Lattice>::isCorner(k) && cidotnormal == 0 ) {

            weightsum += TTraits::Stencil::Weights[idx] * (cidotnormal);

        }
        else {

            distsum += distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] * (cidotnormal);

        }

    }

    static auto model = static_cast<ModelBase<Lattice,TTraits>*>(mModel);

    TParam0thMoment::template get<TTraits::Lattice>(k) =  (mNormalGradient-distsum)/weightsum;

    for (int idx = 1; idx < TTraits::Stencil::Q; idx++) {

        if ( TTraits::Stencil::Ci_x[idx] * normal[0] + TTraits::Stencil::Ci_y[idx] * normal[1] * (TTraits::Lattice::NDIM>1) + TTraits::Stencil::Ci_z[idx] * normal[2] * (TTraits::Lattice::NDIM>2) > 0 ) {
                
            distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = model.computeEquilibrium(k,idx);

        }
        else if ( Geometry<typename TTraits::Lattice>::isCorner(k) && TTraits::Stencil::Ci_x[idx] * normal[0] + TTraits::Stencil::Ci_y[idx] * normal[1] * (TTraits::Lattice::NDIM>1) + TTraits::Stencil::Ci_z[idx] * normal[2] * (TTraits::Lattice::NDIM>2) == 0 ) {
            
            distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = model.computeEquilibrium(k,idx);

        }

    }

    TParam0thMoment::template get<TTraits::Lattice>(k) = oldval;
    
}

template<class TParam0thMoment>
template<class TTraits, class TDistributionType>
inline void Neumann<TParam0thMoment>::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}
