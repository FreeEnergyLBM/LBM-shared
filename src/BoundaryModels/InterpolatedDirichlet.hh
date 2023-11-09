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

        inline void setInterfaceID(int id) {mInterfaceID = id;};

        inline void setInterfaceDistanceFunction(double (*func)(int idx, int k)){

            evalDistanceFunction = func;

        }

    private:

        static double defaultDistanceFunction(int idx, int k) { return true; }

        double (*evalDistanceFunction)(int idx, int k) = &defaultDistanceFunction;

        double mInterfaceVal;
        int mInterfaceID;

};

template<class TTraits, class TDistributionType>
inline void InterpolatedDirichlet::compute(TDistributionType& distribution, int k) {

    if (Geometry<typename TTraits::Lattice>::getBoundaryType(k) != mInterfaceID) return;
    
    for (int idx = 1; idx < TTraits::Stencil::Q; idx++) {

        if (Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == 0 || Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == 6|| Geometry<typename TTraits::Lattice>::getBoundaryType(distribution.streamIndex(k, idx)) == 1) {

            double dist = evalDistanceFunction(distribution.streamIndex(k, idx),distribution.getOpposite(idx));
            
            //std::cout<<dist<<std::endl;
            //#pragma omp critical
            //{
            //std::cout<<dist<<" "<<idx<<" "<<OrderParameter<>::get<typename TTraits::Lattice>(distribution.streamIndex(k, idx))<<std::endl;//fabs(normaldist * normdotci/magci)<<std::endl;
            //}
            //distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = 2.0 * (dist - 1.0) * distribution.getDistributionPointer(k)[distribution.getOpposite(idx)] - (pow(2.0 * dist - 1.0, 2) / (2.0 * dist + 1.0)) * distribution.getDistributionPointer(distribution.streamIndex(k, idx))[distribution.getOpposite(idx)] + 2.0 * ((2.0 * dist - 1.0) / (2.0 * dist + 1.0)) * distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, idx), idx))[idx] + 2.0 * ((3.0 - 2.0 * dist) / (2.0 * dist + 1.0)) * TTraits::Stencil::Weights[idx] * mInterfaceVal;
            //double gamma;
            //distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = 1.0/(1.0+2.0*dist-) * distribution.getDistributionOldPointer(distribution.streamIndex(k, idx))[distribution.getOpposite(idx)] - (pow(2.0 * dist - 1.0, 2) / (2.0 * dist + 1.0)) * distribution.getDistributionOldPointer(distribution.streamIndex(distribution.streamIndex(k, idx), idx))[distribution.getOpposite(idx)] + 2.0 * ((2.0 * dist - 1.0) / (2.0 * dist + 1.0)) * distribution.getDistributionOldPointer(distribution.streamIndex(k, idx))[idx] + 2.0 * ((3.0 - 2.0 * dist) / (2.0 * dist + 1.0)) * TTraits::Stencil::Weights[idx] * mInterfaceVal;
            
            //distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = - distribution.getDistributionOldPointer(distribution.streamIndex(k, idx))[distribution.getOpposite(idx)] + ((2.0 * dist - 1.0) / (2.0 * dist + 1.0)) * distribution.getDistributionOldPointer(distribution.streamIndex(distribution.streamIndex(k, idx), idx))[distribution.getOpposite(idx)] + ((2.0 * dist - 1.0) / (2.0 * dist + 1.0)) * distribution.getDistributionOldPointer(distribution.streamIndex(k, idx))[idx] + 2 * ((2.0) / (2.0 * dist + 1.0)) * TTraits::Stencil::Weights[idx] * mInterfaceVal;

            if (dist <= 0.5) distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -2 * (dist) * distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx)) + (2 * dist - 1) * distribution.getPostCollisionDistribution(distribution.streamIndex(distribution.streamIndex(k, idx), idx),distribution.getOpposite(idx)) + 2 * TTraits::Stencil::Weights[idx] * mInterfaceVal;
            else distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = - 1.0 / (2.0*dist) * distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),distribution.getOpposite(idx)) + (1-1.0 / (2.0 * dist)) * distribution.getPostCollisionDistribution(distribution.streamIndex(k, idx),idx) + 2 * (1.0/(2.0 * dist)) * TTraits::Stencil::Weights[idx] * mInterfaceVal;
            //distribution.getDistributionPointer(distribution.streamIndex(k, idx))[idx] = -distribution.getDistributionOldPointer(distribution.streamIndex(k, idx))[distribution.getOpposite(idx)] + 2*TTraits::Stencil::Weights[idx]*mInterfaceVal;

        }
        
    }    

}

template<class TTraits, class TDistributionType>
inline void InterpolatedDirichlet::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);

}