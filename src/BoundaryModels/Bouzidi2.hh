#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>

class Bouzidi2 : public BoundaryBase {
    public:

        Bouzidi2() : BoundaryBase(5) {}

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);

        inline void setInterfaceDistanceFunction(double (*func)(int idx, int k)){

            evalDistanceFunction = func;

        }

    private:

        static double defaultDistanceFunction(int idx, int k) { return 0.5; }

        double (*evalDistanceFunction)(int idx, int k) = &defaultDistanceFunction;

};

template<class TTraits, class TDistributionType>
inline void Bouzidi2::compute(TDistributionType& distribution, int k) {

    using Lattice = typename TTraits::Lattice;

    if (!this->apply<Lattice>(k)) return;
    
    for (int idx = 1; idx < TTraits::Stencil::Q; idx++) {

        if (this->apply<Lattice>(distribution.streamIndex(k, idx))||Geometry<Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == -1) continue;

        double dist = evalDistanceFunction(distribution.streamIndex(k, idx),distribution.getOpposite(idx));
        
        if (dist <= 0.5) distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = dist * (2 * dist + 1) * distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx)) + (1 + 2 * dist) * (1 - 2 * dist) * distribution.getPostCollisionDistribution(distribution.streamIndex(distribution.streamIndex(k, idx), idx),distribution.getOpposite(idx)) - (dist) * (1 - 2 * dist) * distribution.getPostCollisionDistribution(distribution.streamIndex(distribution.streamIndex(distribution.streamIndex(k, idx), idx), idx),distribution.getOpposite(idx));
        else distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = 1.0 / (dist * (2 * dist + 1)) * distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx)) + (2 * dist - 1) / dist * distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),idx) + (1 - 2 * dist) / (1 + 2 * dist) * distribution.getPostCollisionDistribution(distribution.streamIndex(distribution.streamIndex(k, idx), idx),idx);

}    

}

template<class TTraits, class TDistributionType>
inline void Bouzidi2::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}