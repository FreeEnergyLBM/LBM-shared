#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class DirichletVariable : public BoundaryBase {
    public:

        DirichletVariable() : BoundaryBase(4) {}

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);

        inline void setWallVelocity(const std::unordered_map<int,double>& val) {mInterfaceVal=val;};

    private:

        std::unordered_map<int,double> mInterfaceVal;

};

template<class TTraits, class TDistributionType>
inline void DirichletVariable::compute(TDistributionType& distribution, int k) {

    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

        if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;

        distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx)) + 2*TTraits::Stencil::Weights[idx]*mInterfaceVal[k];

    }    

}

template<class TTraits, class TDistributionType>
inline void DirichletVariable::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}

/*
class DirichletSecondOrder : public BoundaryBase {
    public:

        Dirichlet() { this->setInterfaceID(4); }

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);

        inline void setInterfaceVal(double val) {mInterfaceVal=val;};

    private:

        double mInterfaceVal;

};

template<class TTraits, class TDistributionType>
inline void DirichletSecondOrder::compute(TDistributionType& distribution, int k) {

    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;

    for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {

        if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;

        distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx)) + 2*TTraits::Stencil::Weights[idx]*mInterfaceVal;

    }    

}

template<class TTraits, class TDistributionType>
inline void DirichletSecondOrder::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}
*/