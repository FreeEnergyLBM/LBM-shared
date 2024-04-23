#pragma once
#pragma once
#include "../Parameters.hh"
#include "../Lattice.hh"
#include "../Service.hh"
#include "ForceBase.hh"
#include <math.h>
#include<iostream>
#include<utility>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

template<class TMethod>
class EvaporationPhaseSource : public ForceBase<TMethod> {
    
    public:

        template<class TTraits>
        inline double compute(int k) const;

        template<class TTraits>
        inline double computeXYZ(int xyz, int k) const;

        template<class TTraits>
        inline double computeDensitySource(int k) const;

        inline void setDensityLiquid(double density){

            mDensity = density;

        }

    private:

        double mDensity=1.0;

};

template<class TMethod>
template<class TTraits>
inline double EvaporationPhaseSource<TMethod>::compute(int k) const {

    return -MassSink<>::get<typename TTraits::Lattice>(k) / mDensity;

}

template<class TMethod>
template<class TTraits>
inline double EvaporationPhaseSource<TMethod>::computeXYZ(int xyz, int k) const {

    return 0;//(OrderParameter<>::get<typename TTraits::Lattice>(k) 
            //* Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k,xyz)
            //- OrderParameterOld<>::get<typename TTraits::Lattice>(k) 
            //* VelocityOld<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k,xyz)) / TTraits::Lattice::DT;

}

template<class TMethod>
template<class TTraits>
inline double EvaporationPhaseSource<TMethod>::computeDensitySource(int k) const {

    return - 0.5 * TTraits::Lattice::DT * MassSink<>::get<typename TTraits::Lattice>(k) / mDensity;

}

template<class TMethod>
class EvaporationPressureSource : public ForceBase<TMethod> {
    
    public:

        template<class TTraits>
        inline double compute(int k) const;

        template<class TTraits>
        inline double computeDensitySource(int k) const;

        inline void setLiquidDensity(double density){

            mDensityLiquid = density;
            mPrefactor = 1.0 / (mDensityGas / (1.0 - mInterfaceHumidity)) - 1.0 / mDensityLiquid;

        }

        inline void setGasDensity(double density){

            mDensityGas = density;
            mPrefactor = 1.0 / (mDensityGas / (1.0 - mInterfaceHumidity)) - 1.0 / mDensityLiquid;

        }

        inline void setInterfaceHumidity(double humidity){

            mInterfaceHumidity = humidity;
            mPrefactor = 1.0 / (mDensityGas / (1.0 - mInterfaceHumidity)) - 1.0 / mDensityLiquid;

        }

    private:

        double mPrefactor = 1.0;
        double mDensityLiquid = 1.0;
        double mDensityGas = 1.0;
        double mInterfaceHumidity = 1.0;

};

template<class TMethod>
template<class TTraits>
inline double EvaporationPressureSource<TMethod>::compute(int k) const {

    return TTraits::Stencil::Cs2*Density<>::get<typename TTraits::Lattice>(k)*MassSink<>::get<typename TTraits::Lattice>(k)*mPrefactor;

}

template<class TMethod>
template<class TTraits>
inline double EvaporationPressureSource<TMethod>::computeDensitySource(int k) const {

    return TTraits::Lattice::DT * TTraits::Stencil::Cs2 * 0.5 * (1 - TTraits::Stencil::Weights[0]) * Density<>::get<typename TTraits::Lattice>(k) * MassSink<>::get<typename TTraits::Lattice>(k) * mPrefactor;

}

template<class TMethod>
class EvaporationHumiditySource : public ForceBase<TMethod> {
    
    public:

        template<class TTraits>
        inline double compute(int k) const;

        template<class TTraits>
        inline double computeXYZ(int xyz, int k) const;

        template<class TTraits>
        inline double computeDensitySourceMultiplicative(int k) const;

        inline void setLiquidDensity(double density){

            mDensityLiquid = density;
            mPrefactor = 1.0 / (mDensityGas / (1.0 - mInterfaceHumidity)) - 1.0 / mDensityLiquid;

        }

        inline void setGasDensity(double density){

            mDensityGas = density;
            mPrefactor = 1.0 / (mDensityGas / (1.0 - mInterfaceHumidity)) - 1.0 / mDensityLiquid;

        }

        inline void setInterfaceHumidity(double humidity){

            mInterfaceHumidity = humidity;
            mPrefactor = 1.0 / (mDensityGas / (1.0 - mInterfaceHumidity)) - 1.0 / mDensityLiquid;

        }

    private:

        double mPrefactor = 1.0;
        double mDensityLiquid = 1.0;
        double mDensityGas = 1.0;
        double mInterfaceHumidity = 1.0;

};

template<class TMethod>
template<class TTraits>
inline double EvaporationHumiditySource<TMethod>::compute(int k) const {

    return Humidity<>::get<typename TTraits::Lattice>(k)*MassSink<>::get<typename TTraits::Lattice>(k)*mPrefactor;

}

template<class TMethod>
template<class TTraits>
inline double EvaporationHumiditySource<TMethod>::computeXYZ(int xyz, int k) const {

    return (Humidity<>::get<typename TTraits::Lattice>(k) 
            * Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k,xyz)
            - HumidityOld<>::get<typename TTraits::Lattice>(k) 
            * VelocityOld<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k,xyz)) / TTraits::Lattice::DT;

}

template<class TMethod>
template<class TTraits>
inline double EvaporationHumiditySource<TMethod>::computeDensitySourceMultiplicative(int k) const {

    return 1 - TTraits::Lattice::DT * 0.5 * MassSink<>::get<typename TTraits::Lattice>(k) * mPrefactor;

}