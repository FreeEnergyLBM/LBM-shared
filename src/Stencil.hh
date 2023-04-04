#ifndef STENCIL_HEADER
#define STENCIL_HEADER
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;

//Stencil.hh: This file specifies the stencils that can be used. This determines the discretisation of the model.
//The "D" refers to the number of cartesian dimensions in the stencil, so "D2" would be a 2D simulation. The "Q"
//refers to the number of discrete velocity directions, including the 0 direction. Convention for the directions
//of these discrete velocities is standard in LBM.
//
//Each class contains the information for "D", "Q", the speed of sound squared in the lattice "Cs2", the
//velocity vector arrays "Ci_x[Q]"..., a function to return the x, y or z vector depending on the direction
//passed to it "Ci_xyz(d)", an array of the opposing velocity vector, an array of lattice weights used in the
//LBM collision step "Weights[Q]", the modulus of velocity vectors depending on the direciton "CModulus", a
//currently not implemented array of MRT moments "Moments[Q]" and a function to calculate MRT relaxation times
//"MRTWeights()".
#pragma omp begin declare target
struct D2Q9{ //Most commonly used 2D stencil
    
    static const int D=2; //Number of cartesian directions
    
    static const int Q=9; //Number of velocity directions
    static constexpr double Cs2=0.33333333333333; //Speed of sound squared
    
    static constexpr int Ci_x[Q]={0,1,-1,0,0,1,-1,1,-1}; //Vectors of velocity directions
    static constexpr int Ci_y[Q]={0,0,0,1,-1,1,-1,-1,1}; //There is no convecntion for the ordering of these
    static constexpr int Ci_z[Q]={0,0,0,0,0,0,0,0,0}; //0 array because there is no z direction
    
    enum{x=0,y=1,z=2};
    static auto Ci_xyz(const int d)->const int(&)[Q]{ //Returns velocity direction vector depending on input d
        if (d==x) {
            return Ci_x;
        }
        else if (d==y) {
            return Ci_y;
        }
        else if (d==z) {
            return Ci_z;
        }
        else{
            throw runtime_error(std::string("Error when indexing velocity stencil. Indices must be greater than 0 and less than "+std::to_string(D)));
        }
    }
    static constexpr int Opposites[Q]={0,2,1,4,3,6,5,8,7}; //Opposite vector at a given index
    
    static constexpr double Weights[Q]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0}; //Lattice weights
    
    template<int idx>
    static constexpr int CModulus=Ci_x[idx]*Ci_x[idx]+Ci_y[idx]*Ci_y[idx]; //Returns the modulus of the velocity vector at a given index, used for the MRT weight calculation

    template<int idx> //CURRENTLY NOT IMPLEMENTED BUT WILL GENERATE AN MRT MATRIX FOR EACH STENCIL
    static constexpr int Moments[Q]={1,
                                    3*(CModulus<idx>)-4,
                                    (9*(CModulus<idx>*CModulus<idx>)-21*(CModulus<idx>)+8)/2,
                                    Ci_x[idx],
                                    (3*(CModulus<idx>)-5)*Ci_x[idx],
                                    Ci_y[idx],
                                    (3*(CModulus<idx>)-5)*Ci_y[idx],
                                    Ci_x[idx]*Ci_x[idx]-Ci_y[idx]*Ci_y[idx],
                                    Ci_x[idx]*Ci_y[idx]};
    inline static vector<double> MRTWeights(const double& invtau){ //MRT relaxation rates
        return {0,1,1,0,1,0,1,invtau,invtau};
    }
};
#pragma omp end declare target
struct D3Q19{ //Most commonly used 3D stencil
    static constexpr int D=3;

    static constexpr int Q=19;
    static constexpr double Cs2=0.33333333333333;
    static constexpr int Ci_x[Q]={0,1,-1,0,0,0,0,1,-1,1,-1,0,0,0,0,1,-1,1,-1};
    static constexpr int Ci_y[Q]={0,0,0,1,-1,0,0,1,-1,-1,1,1,-1,1,-1,0,0,0,0};
    static constexpr int Ci_z[Q]={0,0,0,0,0,1,-1,0,0,0,0,1,-1,-1,1,1,-1,-1,1};

    static auto Ci_xyz(const int d)->const int(&)[Q]{
        if (d==0) {
            return Ci_x;
        }
        else if (d==1) {
            return Ci_y;
        }
        else if (d==2) {
            return Ci_z;
        }
        else{
            throw runtime_error(std::string("Error when indexing velocity stencil. Indices must be greater than 0 and less than "+std::to_string(D)));
        }
    }

    static constexpr int Opposites[Q]={0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17};
    static constexpr double Weights[Q]={1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
				                        1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0,
				                        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
				                        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
				                        1.0/36.0, 1.0/36.0, 1.0/36.0};
    template<int idx>
    static const int CModulus=Ci_x[idx]*Ci_x[idx]+Ci_y[idx]*Ci_y[idx]+Ci_z[idx]*Ci_z[idx];
    //static constexpr int CMod[Q]={CMod<0>,CMod<1>,CMod<2>,CMod<3>,CMod<4>,CMod<5>,CMod<6>,CMod<7>,CMod<8>,CMod<9>};

    template<int idx>
    static constexpr int Moments[Q]={1,
                                    19*(CModulus<idx>)-30,
                                    (21*(CModulus<idx>*CModulus<idx>)-53*(CModulus<idx>)+24)/2,
                                    Ci_x[idx],
                                    (5*(CModulus<idx>)-9)*Ci_x[idx],
                                    Ci_y[idx],
                                    (5*(CModulus<idx>)-9)*Ci_y[idx],
                                    Ci_z[idx],
                                    (5*(CModulus<idx>)-9)*Ci_z[idx],
                                    3*Ci_x[idx]*Ci_x[idx]-CModulus<idx>,
                                    (3*CModulus<idx>-5)*(3*Ci_x[idx]*Ci_x[idx]-CModulus<idx>),
                                    Ci_y[idx]*Ci_y[idx]-Ci_z[idx]*Ci_z[idx],
                                    (3*CModulus<idx>-5)*(Ci_y[idx]*Ci_y[idx]-Ci_z[idx]*Ci_z[idx]),
                                    Ci_x[idx]*Ci_y[idx],
                                    Ci_y[idx]*Ci_z[idx],
                                    Ci_x[idx]*Ci_z[idx],
                                    Ci_x[idx]*(Ci_y[idx]*Ci_y[idx]-Ci_z[idx]*Ci_z[idx]),
                                    Ci_y[idx]*(Ci_z[idx]*Ci_z[idx]-Ci_x[idx]*Ci_x[idx]),
                                    Ci_z[idx]*(Ci_x[idx]*Ci_x[idx]-Ci_y[idx]*Ci_y[idx])};
    inline static vector<double> MRTWeights(const double& invtau){
        return {0,1,1,0,1,0,1,0,1,invtau,1,invtau,1,invtau,invtau,invtau,1,1,1};
    }
};

//BELOW: These are just other stencils that I haven't implemented yet

/*
template<typename T>
struct D2Q5{
    static const int D=2;
    static const int Q=5;
    static constexpr double Cs2=0.57735026919;
    static constexpr int Ci_x[Q]={0,1,-1,0,0};
    static constexpr int Ci_y[Q]={0,0,0,1,-1};
    static constexpr int Ci_z[Q]={0,0,0,0,0};
    static auto Ci_xyz(const int d)->const int(&)[Q]{
        if (d==0) {
            return Ci_x;
        }
        if (d==1) {
            return Ci_y;
        }
        if (d==2) {
            return Ci_z;
        }
    }
    static constexpr int Opposites[Q]={0,2,1,4,3};
    static constexpr T Weights[Q]={1/3.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0};
    template<int idx>
    static const int CModulus=Ci_x[idx]*Ci_x[idx]+Ci_y[idx]*Ci_y[idx];
    //static constexpr int CMod[Q]={CMod<0>,CMod<1>,CMod<2>,CMod<3>,CMod<4>,CMod<5>,CMod<6>,CMod<7>,CMod<8>,CMod<9>};

    template<int idx>
    static constexpr int Moments[Q]={1,
                                    3*(CModulus<idx>)-4,
                                    Ci_x[idx],
                                    Ci_y[idx],
                                    Ci_x[idx]*Ci_x[idx]-Ci_y[idx]*Ci_y[idx]};
    inline static vector<int> MRTWeights(const T& invtau){
        return {0,0,0,1,invtau};
    }
};

template<typename T>
struct D3Q19{
    static constexpr int D=3;
    static constexpr int Q=19;
    static constexpr double Cs2=0.57735026919;
    static constexpr int Ci_x[Q]={0,1,-1,0,0,0,0,1,-1,1,-1,0,0,0,0,1,-1,1,-1};
    static constexpr int Ci_y[Q]={0,0,0,1,-1,0,0,1,-1,-1,1,1,-1,1,-1,0,0,0,0};
    static constexpr int Ci_z[Q]={0,0,0,0,0,1,-1,0,0,0,0,1,-1,-1,1,1,-1,-1,1};
    static constexpr int Opposites[Q]={0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17};
    static constexpr T Weights[Q]={1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
				                        1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0,
				                        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
				                        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
				                        1.0/36.0, 1.0/36.0, 1.0/36.0};
    template<int idx>
    static const int CModulus=Ci_x[idx]*Ci_x[idx]+Ci_y[idx]*Ci_y[idx]+Ci_z[idx]*Ci_z[idx];
    //static constexpr int CMod[Q]={CMod<0>,CMod<1>,CMod<2>,CMod<3>,CMod<4>,CMod<5>,CMod<6>,CMod<7>,CMod<8>,CMod<9>};

    template<int idx>
    static constexpr int Moments[Q]={1,
                                    19*(CModulus<idx>)-30,
                                    (21*(CModulus<idx>*CModulus<idx>)-53*(CModulus<idx>)+24)/2,
                                    Ci_x[idx],
                                    (5*(CModulus<idx>)-9)*Ci_x[idx],
                                    Ci_y[idx],
                                    (5*(CModulus<idx>)-9)*Ci_y[idx],
                                    Ci_z[idx],
                                    (5*(CModulus<idx>)-9)*Ci_z[idx],
                                    3*Ci_x[idx]*Ci_x[idx]-CModulus<idx>,
                                    (3*CModulus<idx>-5)*(3*Ci_x[idx]*Ci_x[idx]-CModulus<idx>),
                                    Ci_y[idx]*Ci_y[idx]-Ci_z[idx]*Ci_z[idx],
                                    (3*CModulus<idx>-5)*(Ci_y[idx]*Ci_y[idx]-Ci_z[idx]*Ci_z[idx]),
                                    Ci_x[idx]*Ci_y[idx],
                                    Ci_y[idx]*Ci_z[idx],
                                    Ci_x[idx]*Ci_z[idx],
                                    Ci_x[idx]*(Ci_y[idx]*Ci_y[idx]-Ci_z[idx]*Ci_z[idx]),
                                    Ci_y[idx]*(Ci_z[idx]*Ci_z[idx]-Ci_x[idx]*Ci_x[idx]),
                                    Ci_z[idx]*(Ci_x[idx]*Ci_x[idx]-Ci_y[idx]*Ci_y[idx])};
    inline static vector<int> MRTWeights(const T& invtau){
        return {0,1,1,0,1,0,1,0,1,invtau,1,invtau,1,invtau,invtau,invtau,1,1,1};
    }
};

template<typename T>
struct D3Q7{
    static const int D=3;
    static const int Q=7;
    static constexpr double Cs2=1.0/3.0;
    static constexpr int Ci_x[Q]={0,1,-1,0,0,0,0};
    static constexpr int Ci_y[Q]={0,0,0,1,-1,0,0};
    static constexpr int Ci_z[Q]={0,0,0,0,0,1,-1};
    static constexpr int Opposites[Q]={0,2,1,4,3,6,5};
    static constexpr T Weights[Q]={1.0/4.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0};
    template<int idx>
    static const int CModulus=Ci_x[idx]*Ci_x[idx]+Ci_y[idx]*Ci_y[idx]+Ci_z[idx]*Ci_z[idx];
    //static constexpr int CMod[Q]={CMod<0>,CMod<1>,CMod<2>,CMod<3>,CMod<4>,CMod<5>,CMod<6>,CMod<7>,CMod<8>,CMod<9>};

    template<int idx>
    static constexpr int Moments[Q]={1,
                                    Ci_x[idx],
                                    Ci_y[idx],
                                    Ci_z[idx],
                                    6-7*CModulus<idx>,
                                    3*Ci_x[idx]-(CModulus<idx>*CModulus<idx>),
                                    Ci_y[idx]*Ci_y[idx]-Ci_z[idx]*Ci_z[idx]};
    inline static vector<int> MRTWeights(const T& invtau){
        return {0,0,0,0,1,1,invtau};
    }
};
*/
#endif