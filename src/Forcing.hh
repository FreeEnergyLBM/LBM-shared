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
        if constexpr (std::is_same_v<Stencil,Cartesian>) return traits::Lattice::NDIM; 
        else if constexpr (std::is_same_v<Stencil,AllDirections>) return traits::Stencil::Q;
        else return traits::Lattice::NDIM; 
    } 

    template<class traits>
    inline static constexpr auto ForceStencils() -> const std::array<size_t,sizeof...(Stencils)> { 

        constexpr std::array<size_t,sizeof...(Stencils)> Stencilarray = {StencilToInt<traits,Stencils>()...};
        return Stencilarray;
        
    }

};



struct Guo : ForcingBase<Cartesian> {

    std::vector<double> ma_Force;

    static GuoPrefactor Prefactor;
    
    template<class traits, class force>
    inline void precompute(force& f, int k){
        if (ma_Force.size()<traits::Lattice::NDIM) ma_Force.resize(traits::Lattice::NDIM);
        ma_Force[0]=f.template computeXYZ<traits>(0,k);
        if constexpr (traits::Lattice::NDIM>=2) ma_Force[1]=f.template computeXYZ<traits>(1,k);
        if constexpr (traits::Lattice::NDIM==3) ma_Force[2]=f.template computeXYZ<traits>(2,k);
    }
    
    template<class traits>
    inline double compute(int idx, int k) { //Guo forcing
        
        double prefactor = traits::Stencil::Weights[idx]; //Prefactor for Guo forcing

        double ci_dot_velocity = (traits::Stencil::Ci_x[idx] * Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,0));
        if constexpr (traits::Stencil::D>1) ci_dot_velocity += (traits::Stencil::Ci_y[idx] * Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,1));
        if constexpr (traits::Stencil::D>2) ci_dot_velocity += (traits::Stencil::Ci_z[idx] * Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,2));

        double forceterm = prefactor * (((traits::Stencil::Ci_x[idx] - Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,0)) / traits::Stencil::Cs2
                                            + ci_dot_velocity * traits::Stencil::Ci_x[idx] / (traits::Stencil::Cs2 * traits::Stencil::Cs2)) * ma_Force[0]); //Force Calculation
        if constexpr (traits::Stencil::D>1) forceterm += prefactor * (((traits::Stencil::Ci_y[idx] - Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,1)) / traits::Stencil::Cs2
                                                                        + ci_dot_velocity * traits::Stencil::Ci_y[idx] / (traits::Stencil::Cs2 * traits::Stencil::Cs2)) * ma_Force[1]);
        if constexpr (traits::Stencil::D>2) forceterm += prefactor * (((traits::Stencil::Ci_z[idx] - Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,2)) / traits::Stencil::Cs2
                                                                        + ci_dot_velocity * traits::Stencil::Ci_z[idx] / (traits::Stencil::Cs2 * traits::Stencil::Cs2)) * ma_Force[2]);

        return forceterm;

    }
};

GuoPrefactor Guo::Prefactor;