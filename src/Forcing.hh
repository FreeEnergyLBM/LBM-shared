#pragma once
#include <string>
#include <memory>
#include <utility>
#include <vector>
#include <iostream>
#include <iterator>
#include "Lattice.hh"
#include "Global.hh"
#include "Stencil.hh"
#include "Service.hh"
#include "Collide.hh"

template<class ...Stencils>
struct ForcingBase{

    using mt_Stencils = std::tuple<Stencils...>;

    template<class traits,class Stencil>
    inline static constexpr int StencilToInt(){
        if constexpr (std::is_same_v<Stencil,Cartesian>) return traits::Lattice::m_NDIM; 
        else if constexpr (std::is_same_v<Stencil,AllDirections>) return traits::Stencil::Q;
        else return traits::Lattice::m_NDIM; 
    } 

    template<class traits>
    inline static constexpr auto ForceStencils() -> const std::array<size_t,sizeof...(Stencils)> { 

        constexpr std::array<size_t,sizeof...(Stencils)> Stencilarray = {StencilToInt<traits,Stencils>()...};
        return Stencilarray;
        
    }

};


struct Guo : ForcingBase<Cartesian> {

    std::vector<double> ma_Force;
    
    template<class traits, class force>
    inline void precompute(force& f, int k){
        if (ma_Force.size()<traits::Lattice::m_NDIM) ma_Force.resize(traits::Lattice::m_NDIM);
        ma_Force[0]=f.template computeXYZ<traits>(0,k);
        ma_Force[1]=f.template computeXYZ<traits>(1,k);
        if constexpr (traits::Lattice::m_NDIM==3) ma_Force[2]=f.template computeXYZ<traits>(2,k);
    }
    
    template<class traits>
    inline double compute(int idx, int k) { //Guo forcing
        
        double ci_dot_velocity = 0;
        double forceterm = 0;
        
        double prefactor = traits::Stencil::Weights[idx]; //Prefactor for Guo forcing

        for (int xyz=0;xyz<traits::Stencil::D;xyz++){

            ci_dot_velocity += (traits::Stencil::Ci_xyz(xyz)[idx] * Velocity<>::get<typename traits::Lattice,traits::Lattice::m_NDIM>(k,xyz)); //Dot product of discrete velocity vector
                                                                        //with velocity
        }
        
        for (int xyz = 0; xyz <traits::Stencil::D; xyz++) {

            forceterm += prefactor * (((traits::Stencil::Ci_xyz(xyz)[idx] - Velocity<>::get<typename traits::Lattice,traits::Lattice::m_NDIM>(k,xyz)) / traits::Stencil::Cs2
                                + ci_dot_velocity * traits::Stencil::Ci_xyz(xyz)[idx] / (traits::Stencil::Cs2 * traits::Stencil::Cs2)) * ma_Force[xyz]); //Force
                                                                                                        //Calculation
        }
        
        return forceterm;

    }
};