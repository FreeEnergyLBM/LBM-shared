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

    const static GuoPrefactor Prefactor;
    
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

struct WellBalancedForce : ForcingBase<Cartesian> {

    Guo m_Guo;

    const static GuoPrefactor Prefactor;

    double m_ForceDotVelocity=0;
    
    template<class traits, class force>
    inline void precompute(force& f, int k) {
        m_Guo.precompute<traits>(f,k);
        for (int xyz=0;xyz<traits::Lattice::NDIM;xyz++) m_ForceDotVelocity+=GradientDensity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,xyz)*Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,xyz);
    }
    
    template<class traits>
    inline double compute(int idx, int k) { //Guo forcing
        
        double ci_dot_velocity = 0;
        double ci_dot_ci = 0;
        double forceterm = 0;
        
        double prefactor = traits::Stencil::Weights[idx]; //Prefactor for Guo forcing

        for (int xyz=0;xyz<traits::Stencil::D;xyz++){

            ci_dot_velocity += (traits::Stencil::Ci_xyz(xyz)[idx] *Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,xyz)); //Dot product of discrete velocity vector
                                                                        //with velocity
            ci_dot_ci += traits::Stencil::Ci_xyz(xyz)[idx]*traits::Stencil::Ci_xyz(xyz)[idx];
       
        }

        for (int xyz = 0; xyz <traits::Stencil::D; xyz++) {

            forceterm += GradientDensity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,xyz)*ci_dot_velocity*traits::Stencil::Ci_xyz(xyz)[idx]/traits::Stencil::Cs2; //Force
                                                                                                        //Calculation
        }
        //if(forceterm>0)std::cout<<m_Guo.compute(idx,k)<<" "<<forceterm<<std::endl;
        return m_Guo.compute<traits>(idx,k)+prefactor*forceterm+prefactor * (-m_ForceDotVelocity+0.5*((ci_dot_ci)/traits::Stencil::Cs2-traits::Lattice::NDIM)*m_ForceDotVelocity);

    }
};

struct AllenCahnSourceMethod : ForcingBase<Cartesian> {

    std::vector<double> ma_Force;
    
    const static NoTauDependence Prefactor;

    template<class traits, class force>
    inline void precompute(force& f, int k){
        ma_Force.push_back(f.template computeXYZ<traits>(0,k));
        ma_Force.push_back(f.template computeXYZ<traits>(1,k));
        if constexpr (traits::Lattice::NDIM==3) ma_Force.push_back(f.template computeXYZ<traits>(2,k));
    }
    
    template<class traits>
    inline double compute(int idx, int k) { //Guo forcing
        
        double forceterm = 0;
        
        double prefactor = traits::Stencil::Weights[idx]; //Prefactor for Guo forcing
        
        for (int xyz = 0; xyz < traits::Stencil::D; xyz++) {

            forceterm += prefactor * (traits::Stencil::Ci_xyz(xyz)[idx]) * ma_Force[xyz]; //Force
                                                                                                        //Calculation
        }

        return forceterm;

    }
};

struct He : ForcingBase<Cartesian> {

    std::vector<double> ma_Force;

    const static GuoPrefactor Prefactor;
    
    double velocity_dot_force = 0;

    template<class traits, class force>
    inline void precompute(force& f, int k){
        ma_Force.push_back(f.template computeXYZ<traits>(0,k));
        ma_Force.push_back(f.template computeXYZ<traits>(1,k));
        if constexpr (traits::Lattice::NDIM==3) ma_Force.push_back(f.template computeXYZ<traits>(2,k));

        for (int xyz=0;xyz<traits::Stencil::D;xyz++){

            velocity_dot_force += (ma_Force[xyz] * Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,xyz)); //Dot product of discrete velocity vector
                                                                        //with velocity
        }
    }
    
    template<class traits>
    inline double compute(int idx, int k) { //Guo forcing
        
        double prefactor = CollisionBase<typename traits::Lattice,typename traits::Stencil>::computeGamma(&Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,0),idx); //Prefactor for Guo forcing

        double ci_dot_force=0;

        for (int xyz=0;xyz<traits::Stencil::D;xyz++){

            ci_dot_force += (traits::Stencil::Ci_xyz(xyz)[idx] * ma_Force[xyz]); //Dot product of discrete velocity vector
                                                                        //with velocity
        }

        return prefactor*(ci_dot_force-velocity_dot_force);

    }
};


struct Lee : ForcingBase<Cartesian,AllDirections> {

    std::vector<double> ma_Force;
    std::vector<double> ma_ForceQ;

    const static NoTauDependence Prefactor;
    
    double velocity_dot_force = 0;

    template<class traits, class force>
    inline void precompute(force& f, int k){
        ma_Force.push_back(f.template computeXYZ<traits>(0,k));
        ma_Force.push_back(f.template computeXYZ<traits>(1,k));
        if constexpr (traits::Lattice::NDIM==3) ma_Force.push_back(f.template computeXYZ<traits>(2,k));
        for (int q = 0; q < traits::Stencil::Q; q++) ma_ForceQ.push_back(f.template computeQ<traits>(q,k));
        
        for (int xyz=0;xyz<traits::Stencil::D;xyz++){

            velocity_dot_force += (ma_Force[xyz] * Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,xyz)); //Dot product of discrete velocity vector
                                                                        //with velocity
        }
        
    }
    
    template<class traits>
    inline double compute(int idx, int k) { //Guo forcing
        
        double prefactor = traits::Lattice::DT * CollisionBase<typename traits::Lattice,typename traits::Stencil>::computeGamma(&Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,0),idx); //Prefactor for Guo forcing

        return prefactor*(ma_ForceQ[idx]-velocity_dot_force);

    }
};


struct LeeGamma0 : ForcingBase<Cartesian,AllDirections> {

    std::vector<double> ma_Force;
    std::vector<double> ma_ForceQ;

    const static NoTauDependence Prefactor;
    
    double velocity_dot_force = 0;

    template<class traits, class force>
    inline void precompute(force& f, int k){
        ma_Force.push_back(f.template computeXYZ<traits>(0,k));
        ma_Force.push_back(f.template computeXYZ<traits>(1,k));
        if constexpr (traits::Lattice::NDIM==3) ma_Force.push_back(f.template computeXYZ<traits>(2,k));
        for (int q = 0; q < traits::Stencil::Q; q++) ma_ForceQ.push_back(f.template computeQ<traits>(q,k));

        for (int xyz=0;xyz<traits::Stencil::D;xyz++){

            velocity_dot_force += (ma_Force[xyz] * Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,xyz)); //Dot product of discrete velocity vector
                                                                        //with velocity
        }
    }
    
    template<class traits>
    inline double compute(int idx, int k) { //Guo forcing
        
        double prefactor = traits::Lattice::DT * traits::Stencil::Weights[idx]; //Prefactor for Guo forcing

        return prefactor*(ma_ForceQ[idx]-velocity_dot_force);

    }
};


struct LeeMuLocal : ForcingBase<AllDirections> {

    std::vector<double> ma_ForceQ;

    const static NoTauDependence Prefactor;

    template<class traits, class force>
    const inline void precompute(force& f, int k){

        for (int q = 0; q < traits::Stencil::Q; q++) ma_ForceQ.push_back(f.template computeQ<traits>(q,k));

    }

    template<class traits, class force>
    inline void precompute(force& f, int q, int k){

        ma_ForceQ.push_back(f.template computeQ<traits>(q,k));

    }
    
    template<class traits>
    inline double compute(int idx, int k) { //Guo forcing

        double prefactor = 0.5 * traits::Lattice::DT * CollisionBase<typename traits::Lattice,typename traits::Stencil>::computeGamma(&Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(k,0),idx); //Prefactor for Guo forcing

        return prefactor*(ma_ForceQ[idx]);

    }
};


struct LeeMuNonLocal : ForcingBase<AllDirections> {

    std::vector<double> ma_ForceQ;

    const static NoTauDependence Prefactor;

    template<class traits, class force>
    const inline void precompute(force& f, int k){

        for (int q = 0; q < traits::Stencil::Q; q++) ma_ForceQ.push_back(f.template computeQ<traits>(q,k));

    }

    template<class traits, class force>
    const inline void precompute(force& f, int q, int k){

        ma_ForceQ.push_back(f.template computeQ<traits>(q,k));

    }
    
    template<class traits>
    const inline double compute(int idx, int k) { //Guo forcing

        using data = Data_Base<typename traits::Lattice, typename traits::Stencil>;
        
        const double prefactor = 0.5 * traits::Lattice::DT * CollisionBase<typename traits::Lattice,typename traits::Stencil>::computeGamma(&Velocity<>::get<typename traits::Lattice,traits::Lattice::NDIM>(data::getInstance().getNeighbors()[k * traits::Stencil::Q+idx],0),idx); //Prefactor for Guo forcing

        return prefactor*(ma_ForceQ[idx]);

    }
};