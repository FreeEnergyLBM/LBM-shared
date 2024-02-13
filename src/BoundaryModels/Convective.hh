#pragma once
#include "../Parameters.hh"
#include "BoundaryBase.hh"
#include<iostream>


class Convective : public BoundaryBase {
    public:

        Convective() { this->setNodeID(4, true); } // TMP: Default NodeID warning

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);


};

template<class TTraits, class TDistributionType>
inline void Convective::compute(TDistributionType& distribution, int k) { //CHANGE THIS SO YOU DONT NEED TO COMMUNICATE
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    const int& normalq = Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;
    const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;
    //#pragma omp critical
    //{
    //std::cout<<normalq<<std::endl;
    //}
    //std::cout<<normalq<<std::endl;

    for (int idx = 1; idx < Stencil::Q; idx++) {
        
        //if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;
        if (BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(distribution.streamIndex(k, idx)).Id!=0) continue;

        double normalvelocity = 0;
        double magnormal = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            normalvelocity += -normal[xyz]*Velocity<>::get<Lattice, Lattice::NDIM>(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq),xyz);
            magnormal += pow(normal[xyz],2);
        }

        magnormal = sqrt(magnormal);

        normalvelocity *= 1./magnormal;

        distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = (distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[idx]+normalvelocity*distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[idx])/(1+normalvelocity);
        //#pragma omp critical
        //{
        //std::cout<<k<<" "<<(distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[idx])<<std::endl;
        //}
        

    }    

}

template<class TTraits, class TDistributionType>
inline void Convective::communicate(TDistributionType& distribution) {

    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
    Lattice::communicateDistributionAllOld(distribution);

}

template<class TForceTuple>
class Convective2 : public BoundaryBase {
    public:

        Convective2() { this->setNodeID(4, true); } // TMP: Default NodeID warning

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        inline void setForceTuple(const TForceTuple& tup) {mt_Forces = tup;}

        template<class TTraits, class TDistributionType>
        inline void communicate(TDistributionType& mDistribution);

        template<class TTraits>
        inline void communicateProcessor();

        template<class TTraits>
        inline void runProcessor(int k);

        inline void setVelocityCalculator(double (*v)(const double* distribution, TForceTuple& forcetuple,  const double& density, int xyz, int k)) {mVelocityCalculator=v;}

    private:

        TForceTuple mt_Forces;

        static double defaultVelocityCalculator(const double* distribution, TForceTuple& forcetuple, const double& density, int xyz, int k) { return 0; }

        double (*mVelocityCalculator)(const double* distribution, TForceTuple& forcetuple, const double& density, int xyz, int k) = &defaultVelocityCalculator;

        double mVelocity = 0;
        int mCount = 0;


};

template<class TForceTuple>
template<class TTraits>
inline void Convective2<TForceTuple>::runProcessor(int k) { //CHANGE THIS SO YOU DONT NEED TO COMMUNICATE
    
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    using DataType = Data_Base<typename TTraits::Lattice, typename TTraits::Stencil>;

    DataType& data = DataType::getInstance();

    const int& normalq = Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;
    const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;

    double normalvelocity=0;
    double magnormal = 0;

    for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
        normalvelocity += normal[xyz]*(Velocity<>::get<Lattice, Lattice::NDIM>(data.getNeighbor(data.getNeighbor(k,normalq),normalq),xyz));//Velocity<>::get<Lattice, Lattice::NDIM>(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq),xyz);
        magnormal += pow(normal[xyz],2);
    }

    magnormal = sqrt(magnormal);

    normalvelocity *= 1./magnormal;
    /*
    #pragma omp parallel reduction (+:mVelocity,mCount)
    {
        mVelocity+=normalvelocity;
        mCount+=1;
    }
    */
    #pragma omp critical
    {
        mVelocity+=normalvelocity;
        mCount+=1;
    }
    //std::cout<<"HERE "<<mVelocity/((double)mCount)<<std::endl;
    //#pragma omp critic
    //if (fabs(normalvelocity)>fabs(mVelocity)) mVelocity=-normalvelocity;
    
}

template<class TForceTuple>
template<class TTraits, class TDistributionType>
inline void Convective2<TForceTuple>::compute(TDistributionType& distribution, int k) { //CHANGE THIS SO YOU DONT NEED TO COMMUNICATE
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    if (!this->apply<Lattice>(k)) return;

    const int& normalq = Stencil::QMap.find(BoundaryLabels<TTraits::Lattice::NDIM>::template get<Lattice>(k).NormalDirection)->second;
    const std::array<int8_t,TTraits::Lattice::NDIM>& normal = BoundaryLabels<TTraits::Lattice::NDIM>::template get<typename TTraits::Lattice>(k).NormalDirection;
    
    for (int idx = 1; idx < Stencil::Q; idx++) {
        
        //if (this->apply<Lattice>(distribution.streamIndex(k, idx))) continue;
        double normdotci=0;
        for (int xyz = 0; xyz<Lattice::NDIM; xyz++){
            normdotci +=  normal[xyz] * Stencil::Ci_xyz(xyz)[idx];
        }
        if (normdotci<=0) continue;
        //std::cout<<idx<<std::endl;
        
        double civelocity = 0;
        double normalvelocity = 0;
        double magnormal = 0;
        for (int xyz = 0; xyz < TTraits::Lattice::NDIM; xyz++) {
            civelocity += 3*Stencil::Ci_xyz(xyz)[idx]*mVelocityCalculator(distribution.getDistributionPointer(distribution.streamIndex(k, normalq)),mt_Forces,Density<>::get<Lattice>(distribution.streamIndex(k, normalq)),xyz,distribution.streamIndex(k, normalq))-4*Stencil::Ci_xyz(xyz)[idx]*mVelocityCalculator(distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq)),mt_Forces,Density<>::get<Lattice>(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq)),xyz,distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))+Stencil::Ci_xyz(xyz)[idx]*mVelocityCalculator(distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq), normalq)),mt_Forces,Density<>::get<Lattice>(distribution.streamIndex(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq), normalq)),xyz,distribution.streamIndex(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq), normalq));//Velocity<>::get<Lattice, Lattice::NDIM>(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq),xyz);
            magnormal += pow(normal[xyz],2);
            normalvelocity += normal[xyz]*mVelocityCalculator(distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq)),mt_Forces,Density<>::get<Lattice>(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq)),xyz,distribution.streamIndex(distribution.streamIndex(k, normalq), normalq));//Velocity<>::get<Lattice, Lattice::NDIM>(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq),xyz);
            
        }

        magnormal = sqrt(magnormal);

        normalvelocity *= 1./magnormal;
        //if (idx==6) distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[8] = (distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[8]+normalvelocity*distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[8])/(1+normalvelocity);
        //else if (idx==8) distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[6] = (distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[6]+normalvelocity*distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[6])/(1+normalvelocity);
        //distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = (distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[idx]+normalvelocity*distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[idx])/(1+normalvelocity);
        
        //distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = (distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[idx]+mVelocity*distribution.getDistributionPointer(distribution.streamIndex(distribution.streamIndex(k, normalq), normalq))[idx])/(1+mVelocity);
        //distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[idx]+3*TTraits::Stencil::Weights[idx]*mVelocity/(double)mCount/2.0*civelocity;
        distribution.getDistributionPointer(distribution.streamIndex(k, normalq))[idx] = distribution.getDistributionOldPointer(distribution.streamIndex(k, normalq))[idx]+3*TTraits::Stencil::Weights[idx]*normalvelocity/2.0*civelocity;

        //#pragma omp critical
        //{
        //if (TIME==50000) std::cout<<mCount<<std::endl;
        //}

    }  

}

template<class TForceTuple>
template<class TTraits, class TDistributionType>
inline void Convective2<TForceTuple>::communicate(TDistributionType& distribution) {
    
    using Lattice = typename TTraits::Lattice;
    Lattice::communicateDistributionAll(distribution);
    Lattice::communicateDistributionAllOld(distribution);

}

template<class TForceTuple>
template<class TTraits>
inline void Convective2<TForceTuple>::communicateProcessor() {
    
    #pragma omp single
    {
        //if(TIME%10000==0)std::cout<<"HERE "<<mVelocity<<std::endl;
    
    
    mVelocity=0;
    mCount=0;
    }
    
}
