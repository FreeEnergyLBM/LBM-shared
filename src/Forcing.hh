#pragma once
#include <string>
#include <memory>
#include <utility>
#include <vector>
#include <iostream>
#include "Lattice.hh"
#include "Global.hh"
#include "Stencil.hh"

template<class lattice, class stencil>
struct Guo{
    Guo(){
        ma_Force[0]=0.;
        ma_Force[1]=0.;
        if constexpr (lattice::m_NDIM==3) ma_Force[2]=0.;
    }
    Guo(Guo& other) { std::copy(other.ma_Force,other.ma_Force+4,ma_Force); }
    Guo(const Guo& other) { std::copy(other.ma_Force,other.ma_Force+4,ma_Force); }

    InverseTau<lattice> m_InvTau;
    Velocity<lattice> m_Velocity;

    double ma_Force[lattice::m_NDIM];
    
    template<class force>
    const inline void precompute(force& f, int k){
        ma_Force[0]+=f.computeXYZ(0,k);
        ma_Force[1]+=f.computeXYZ(1,k);
        if constexpr (lattice::m_NDIM==3) ma_Force[2]+=f.computeXYZ(2,k);
    }
    const inline double compute(int idx, int k) const { //Guo forcing
        
        double ci_dot_velocity = 0;
        double forceterm = 0;
        
        double prefactor = (1 - lattice::m_DT * m_InvTau.getParameter(k) / 2.0) * stencil::Weights[idx]; //Prefactor for Guo forcing

        for (int xyz=0;xyz<stencil::D;xyz++){

            ci_dot_velocity += (stencil::Ci_xyz(xyz)[idx] * m_Velocity.getParameter()[k * stencil::D]); //Dot product of discrete velocity vector
                                                                        //with velocity
        }
        
        for (int xyz = 0; xyz <stencil::D; xyz++) {

            forceterm += prefactor * (((stencil::Ci_xyz(xyz)[idx] - m_Velocity.getParameter()[k * stencil::D]) / stencil::Cs2
                                + ci_dot_velocity * stencil::Ci_xyz(xyz)[idx] / (stencil::Cs2 * stencil::Cs2)) * ma_Force[xyz]); //Force
                                                                                                        //Calculation
        }

        return forceterm;

    }
};

struct ForcingNone{
    template<class force>
    const inline void precompute(force& f, int k) {}
    const inline double compute(int idx, int k) { return 0.; }
};
