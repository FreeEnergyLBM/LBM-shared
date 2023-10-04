#pragma once
#include "../Service.hh"
#include "GradientBase.hh"

struct CentralXYZInterfaceNoSolid : InterfaceGradient<Cartesian> {
    
    template<class TTraits, class TParameter>
    inline double compute( int direction, int k, int num = 0);

    template<class TObj>
    using GradientType = Gradient<TObj,TObj::instances>;

    inline void setInterfaceCondition(bool (*condition)(int k)){

        evalInterfaceCondition=condition;

    }

    inline void setInterfaceDistance(double (*distance)(int k, int idx)){

        evalInterfaceDistance = distance;

    }

    inline void setInterfaceVal(double value){

        mInterfaceVal = value;

    }

    double mInterfaceVal = 0;

    static bool defaultCondition(int k) { return true; }

    bool (*evalInterfaceCondition)(int k) = &defaultCondition;

    static double defaultDistance(int k, int idx) { return 0.5; }

    double (*evalInterfaceDistance)(int k, int idx) = &defaultDistance;
    
};

template<class TTraits, class TParameter>
inline double CentralXYZInterfaceNoSolid::compute(int direction, int k, int num) {
    
    using Lattice = typename TTraits::Lattice;
    using Stencil = typename TTraits::Stencil;

    using DataType = Data_Base<Lattice, Stencil>;

    DataType& data = DataType::getInstance();

    double gradientsum = 0;

    if (Geometry<Lattice>::getBoundaryType(k)==0){

        for (int idx = 1; idx <Stencil::Q; idx++) {

            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k,idx))==1)) {
                
                if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==0||Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==4) {
                
                    gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(k, num));
                    //std::cout<<"Minusinterface "<<Stencil::Ci_xyz(direction)[idx]<<" "<<(TParameter::template get<Lattice>(k, num))<<" "<<(TParameter::template get<Lattice>(data.getNeighbor(k, Stencil::Opposites[idx]), num))<<" "<<(TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), Stencil::Opposites[idx]), num))<<std::endl;
                }
                else {

                    double interfacedistance = evalInterfaceDistance(k,Stencil::Opposites[idx]);
                    gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(k, num)) / (1 + interfacedistance);

                }

            }
            else {

                if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))==0 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))==6 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))==4) && (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==0 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==6 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==4)) {

                    gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(data.getNeighbor(k, idx), num));
                    //std::cout<<"Nointerface "<<Stencil::Ci_xyz(direction)[idx]<<" "<<TParameter::template get<Lattice>(data.getNeighbor(k, idx), num)<<std::endl;
                }
                else if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == 5 && (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==0 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==6)) {

                    
                    double interfacedistance = evalInterfaceDistance(k,idx);
                    gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (mInterfaceVal) / (1 + interfacedistance);
                    //std::cout<<"Plusinterface "<<Stencil::Ci_xyz(direction)[idx]<<" "<<evalInterfaceDistance(k,idx)<<" "<<(mInterfaceVal)<<std::endl;

                }
                else if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) != 5 && Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==5){

                    double interfacedistance = evalInterfaceDistance(k,Stencil::Opposites[idx]);
                    //std::cout<<TParameter::template get<Lattice>(data.getNeighbor(k, idx), num)<<" "<<interfacedistance<<std::endl;
                    gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(data.getNeighbor(k, idx), num)) / (1 + interfacedistance);
                    //std::cout<<"Minusinterface "<<Stencil::Ci_xyz(direction)[idx]<<" "<<interfacedistance<<" "<<(TParameter::template get<Lattice>(k, num))<<" "<<(TParameter::template get<Lattice>(data.getNeighbor(k, idx), num))<<" "<<(TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx), num))<<std::endl;
                }                

            }
            
        }

    }
    else if (Geometry<Lattice>::getBoundaryType(k)==4){

        for (int idx = 1; idx <Stencil::Q; idx++) {

            if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k,idx))==1)) {
                
                if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==0||Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==4) {
                
                    gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (0);
                    //std::cout<<"Minusinterface "<<Stencil::Ci_xyz(direction)[idx]<<" "<<(TParameter::template get<Lattice>(k, num))<<" "<<(TParameter::template get<Lattice>(data.getNeighbor(k, Stencil::Opposites[idx]), num))<<" "<<(TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, Stencil::Opposites[idx]), Stencil::Opposites[idx]), num))<<std::endl;
                }
                else {

                    double interfacedistance = evalInterfaceDistance(k,Stencil::Opposites[idx]);
                    gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(k, num)) / (1 + interfacedistance);

                }

            }
            else if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==1)) {
                
                if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))==0 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))==6 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))==4)) {

                    gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(data.getNeighbor(k, idx), num))*0.75;
                    //std::cout<<"Nointerface "<<Stencil::Ci_xyz(direction)[idx]<<" "<<TParameter::template get<Lattice>(data.getNeighbor(k, idx), num)<<std::endl;
                }

            } 
            else {

                if ((Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))==0 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))==6 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx))==4) && (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==0 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==6 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==4)) {

                    gradientsum += Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(data.getNeighbor(k, idx), num));
                    //std::cout<<"Nointerface "<<Stencil::Ci_xyz(direction)[idx]<<" "<<TParameter::template get<Lattice>(data.getNeighbor(k, idx), num)<<std::endl;
                }
                else if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) == 5 && (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==0 || Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==6)) {

                    
                    double interfacedistance = evalInterfaceDistance(k,idx);
                    gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (mInterfaceVal) / (1 + interfacedistance);
                    //std::cout<<"Plusinterface "<<Stencil::Ci_xyz(direction)[idx]<<" "<<evalInterfaceDistance(k,idx)<<" "<<(mInterfaceVal)<<std::endl;

                }
                else if (Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, idx)) != 5 && Geometry<Lattice>::getBoundaryType(data.getNeighbor(k, Stencil::Opposites[idx]))==5){

                    double interfacedistance = evalInterfaceDistance(k,Stencil::Opposites[idx]);
                    //std::cout<<TParameter::template get<Lattice>(data.getNeighbor(k, idx), num)<<" "<<interfacedistance<<std::endl;
                    gradientsum += 2 * Stencil::Weights[idx] * Stencil::Ci_xyz(direction)[idx] * (TParameter::template get<Lattice>(data.getNeighbor(k, idx), num)) / (1 + interfacedistance);
                    //std::cout<<"Minusinterface "<<Stencil::Ci_xyz(direction)[idx]<<" "<<interfacedistance<<" "<<(TParameter::template get<Lattice>(k, num))<<" "<<(TParameter::template get<Lattice>(data.getNeighbor(k, idx), num))<<" "<<(TParameter::template get<Lattice>(data.getNeighbor(data.getNeighbor(k, idx), idx), num))<<std::endl;
                }                

            }
            
        }

    }

    return 1.0 / (Stencil::Cs2 * Lattice::DT) * gradientsum;

}
