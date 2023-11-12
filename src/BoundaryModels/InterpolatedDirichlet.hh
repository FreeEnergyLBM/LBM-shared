#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>

class InterpolatedDirichlet : public BoundaryBase {
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);

        inline void setInterfaceVal(double val) {mInterfaceVal = val;};

            inline void setInterfaceID(int id) {mInterfaceID[0]=id;};

        inline void setInterfaceID(const std::vector<int>& id) {mInterfaceID=id;};

        inline void setInterfaceDistanceFunction(double (*func)(int idx, int k)){

            evalDistanceFunction = func;

        }

    private:

        static double defaultDistanceFunction(int idx, int k) { return true; }

        double (*evalDistanceFunction)(int idx, int k) = &defaultDistanceFunction;

        double mInterfaceVal;
        std::vector<int> mInterfaceID = {4};

};

template<class TTraits, class TDistributionType>
inline void InterpolatedDirichlet::compute(TDistributionType& distribution, int k) {

    for (int i : mInterfaceID){
        if(Geometry<typename TTraits::Lattice>::getBoundaryType(k) == i) goto runloop;
    }

    return;

    runloop:
    
        for (int idx = 1; idx < TTraits::Stencil::Q; idx++) {

            bool cont = true;

            for (int i : mInterfaceID){
                if(Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == i) goto dontapply;
            }

            cont = false;

            dontapply:
                if (cont) continue;

            double dist = evalDistanceFunction(distribution.streamIndex(k, idx),distribution.getOpposite(idx));
            
            if (dist <= 0.5) distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -2 * (dist) * distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx)) + (2 * dist - 1) * distribution.getPostCollisionDistribution(distribution.streamIndex(distribution.streamIndex(k, idx), idx),distribution.getOpposite(idx)) + 2 * TTraits::Stencil::Weights[idx] * mInterfaceVal;
            else distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = - 1.0 / (2.0*dist) * distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx)) + (1-1.0 / (2.0 * dist)) * distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),idx) + 2 * (1.0/(2.0 * dist)) * TTraits::Stencil::Weights[idx] * mInterfaceVal;

    }    

}

template<class TTraits, class TDistributionType>
inline void InterpolatedDirichlet::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}