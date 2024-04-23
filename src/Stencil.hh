#pragma once
#include <stdexcept>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>

/**
 * \file Stencil.hh
 * \brief  This file specifies the stencils that can be used. This determines the velocity discretisation of the model.
 * The "D" refers to the number of cartesian dimensions in the stencil, so "D2" would be a 2D simulation. The "Q"
 * refers to the number of discrete velocity directions, including the 0 direction. Convention for the directions
 * of these discrete velocities is standard in LBM.
 *
 * Each class contains the information for "D", "Q", the speed of sound squared in the lattice "Cs2", the
 * velocity vector arrays "Ci_x[Q]"..., a function to return the x, y or z vector depending on the direction
 * passed to it "Ci_xyz(d)", an array of the opposing velocity vector, an array of lattice weights used in the
 * LBM collision step "Weights[Q]", the modulus of velocity vectors depending on the direciton "CModulus", a
 * currently not implemented array of MRT moments "Moments[Q]" and a function to calculate MRT relaxation times
 * "MRTWeights()".
 */

//using boost::hash_combine
template <class T>
inline void hash_combine(std::size_t& seed, T const& v)
{
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
    template<typename T>
    struct hash<vector<T>>
    {
        typedef vector<T> argument_type;
        typedef std::size_t result_type;
        result_type operator()(argument_type const& in) const
        {
            size_t size = in.size();
            size_t seed = 0;
            for (size_t i = 0; i < size; i++)
                //Combine the hash of the current vector with the hashes of the previous ones
                hash_combine(seed, in[i]);
            return seed;
        }
    };
}

struct Stencil3D{};

struct StencilQ{};

struct StencilBase {};

/**
 * \brief Class to set one MRT matrix row based on the given moment class.
 * \tparam N Element of the MRT matrix row to start from. Also the velocity index.
 * \tparam moment Class containing the moment information. Accepts N as a template paramter.
 * This class will count down from N until 0, setting these elements of the matrix row.
 */
template<int N,template<int> class moment>
struct MRTRowSetter{
    
    /**
     * \brief This function sets one element N of an mrt matrix row idx. It then calls itself for the next element.
     * \param matrix Pointer to the first element of the MRT matrix row idx.
     * \param idx Index of the moment/row to be set.
     */
    static constexpr void setMRTRow(int* matrix,int idx){
        matrix[N]=moment<N>::ma_Moments[idx];
        MRTRowSetter<N-1,moment>::setMRTRow(matrix,idx);
    }

};


/**
 * \brief Class to set the first element in an  MRT matrix row based on the given moment class. Specialisation of MRTRowSetter.
 * \tparam moment Class containing the moment information. Accepts N as a template paramter, although N is 0 in this case.
 */
template<template<int> class moment>
struct MRTRowSetter<0,moment>{
    
    /**
     * \brief This function sets the first element of an mrt matrix row idx.
     * \param matrix Pointer to the first element of the MRT matrix row idx.
     * \param idx Index of the moment/row to be set.
     * This class should be used to terminate the recursive setting of elements in the MRT matrix, counting from N to 0 in this class.
     */
    static constexpr void setMRTRow(int* matrix,int idx){
        matrix[0]=moment<0>::ma_Moments[idx];
    }

};

/**
 * \brief Class to set the MRT matrix for a given stencil.
 * \tparam N Index of the current moment being set.
 * \tparam Q Size of the matrix QxQ
 * \tparam moment Class containing the moment information. Accepts the velocity index as a template paramter.
 * This class will count down from N until 0, setting these rows of the MRT matrix.
 */
template<int N, int Q,template<int> class moment>
struct MRTMatrixSetter{
    
    /**
     * \brief This function sets one row N of an mrt matrix. It then calls itself for the next row, decrimenting N.
     * \param matrix Pointer to the first element of the MRT matrix
     */
    static constexpr void setMRTMatrix(int* matrix){
        MRTRowSetter<Q-1,moment>::setMRTRow(&matrix[Q*N],N);
        MRTMatrixSetter<N-1,Q,moment>::setMRTMatrix(matrix);
    }

};


/**
 * \brief Class to set the first row of an matrix for a given stencil. Specialisation of MRTMatrix Setter.
 * \tparam Q Size of the matrix QxQ.
 * \tparam moment Class containing the moment information. Accepts the velocity index as a template paramter.
 */
template<int Q,template<int> class moment>
struct MRTMatrixSetter<0,Q,moment>{
    
    /**
     * \brief This function sets the first row of an mrt matrix.
     * \param matrix Pointer to the first element of the MRT matrix
     * This is used to terminate the recursive MRTMatrixSetter function.
     */
    static constexpr void setMRTMatrix(int* matrix){
        MRTRowSetter<Q-1,moment>::setMRTRow(&matrix[0],0);
    }

};

/**
 * \brief This function computes and returns the MRT matrix for a stencil.
 * \tparam Q The number of velociy directions in the stencil.
 * \tparam moment Class containing the moment array. Accepts the velocity index as a template argument.
 * \return The full MRT matrix for the stencil.
 */
template<int Q,template<int> class moment>
static constexpr std::array<int,Q*Q> GenerateMRTMatrix(){
    std::array<int,Q*Q> matrix{};
    MRTMatrixSetter<Q-1,Q,moment>::setMRTMatrix(&matrix[0]);
    return matrix;
}

/**
 * \brief 1D velocity stencil class.
 */
struct D1Q3 : StencilBase {
    
    static constexpr int D = 1; //!<Number of cartesian directions.
    static constexpr int Q = 3; //!<Number of velocity directions.
    static constexpr double Cs2 = 1.0/3.0; //!<Speed of sound squared.
    
    static constexpr int Ci_x[Q] = {0, 1, -1}; //!<X component vector of velocity directions.
    static constexpr int Ci_y[Q] = {0, 0, 0}; //!<Y component vector of velocity directions. 0 for 1D stencil.
    static constexpr int Ci_z[Q] = {0, 0, 0}; //!<Z component vector of velocity directions. 0 for 1D stencil.

    static std::map<std::array<int8_t,D>,int> QMap; //!<Map of cartesian vector to velocity index.

    enum{x = 0, y = 1, z = 2};

    /**
     * \brief Returns velocity direction vector depending on input d.
     * \param d Cartesian direction (x=0, y=1, z=2).
     * \return Vector of d component of velocity vectors.
     */
    inline static auto Ci_xyz(const int d) -> const int(&)[Q] {

        if (d == x) {
            return Ci_x;
        }
        else if (d == y) {
            return Ci_y;
        }
        return Ci_z;

    }

    static constexpr int Opposites[Q] = {0, 2, 1}; //!<Opposite vector at a given index.
    
    static constexpr double Weights[Q] = {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}; //!<Lattice weights.

    template<int idx>
    static constexpr int CModulus = Ci_x[idx] * Ci_x[idx]; //!<Computes the modulus of the velocity vector at a given index, used for the MRT weight calculation.
    
    /**
     * \brief Class containing an array of the MRT matrix moments for a given idx.
     * \tparam idx Velocity vector index.
     */
    template<int idx>
    struct Moments{
        static constexpr int ma_Moments[Q]={1,
                                        Ci_x[idx],
                                        3 * (CModulus<idx>)-2};
    };

    static constexpr std::array<int,Q*Q> MRTMatrix = GenerateMRTMatrix<Q,Moments>(); //!<MRT matrix array.

    /**
     * \brief Function that returns the diagonal entries of the MRT weight matrix for a given inverse relaxation time.
     * \param invtau Inverse of the relaxation time.
     */
    inline static std::vector<double> MRTWeights(const double& invtau){
        return {0,0,invtau};
    }

};

std::map<std::array<int8_t,D1Q3::D>,int> D1Q3::QMap = {{{0},0},
                                                      {{1},1},
                                                      {{-1},2}};

/**
 * \brief Simple 2D stencil with no diagonal velocity directions.
 */
struct D2Q5 : StencilBase {
    
    static constexpr int D = 2; //!<Number of cartesian directions.
    static constexpr int Q = 5; //!<Number of velocity directions.
    static constexpr double Cs2 = 1.0/3.0; //Speed of sound squared.
    
    static constexpr int Ci_x[Q] = {0, 1, -1, 0, 0}; //!<X component vector of velocity directions.
    static constexpr int Ci_y[Q] = {0, 0, 0, 1, -1}; //!<Y component vector of velocity directions.
    static constexpr int Ci_z[Q] = {0, 0, 0, 0, 0}; //!<Z component vector of velocity directions. 0 for 2D stencil.

    static std::map<std::array<int8_t,D>,int> QMap; //!<Map of cartesian vector to velocity index.

    enum{x = 0, y = 1, z = 2};
    
    /**
     * \brief Returns velocity direction vector depending on input d.
     * \param d Cartesian direction (x=0, y=1, z=2).
     * \return Vector of d component of velocity vectors.
     */
    inline static auto Ci_xyz(const int d) -> const int(&)[Q] {

        if (d == x) {
            return Ci_x;
        }
        else if (d == y) {
            return Ci_y;
        }
        return Ci_z;

    }

    static constexpr int Opposites[Q] = {0, 2, 1, 4, 3}; //!<Opposite vector at a given index
    
    static constexpr double Weights[Q] = {1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0}; //!<Lattice weights

    template<int idx>
    static constexpr int CModulus = Ci_x[idx] * Ci_x[idx] + Ci_y[idx] * Ci_y[idx]; //!<Computes the modulus of the velocity vector at a given index, used for the MRT weight calculation
    
    /**
     * \brief Class containing an array of the MRT matrix moments for a given idx.
     * \tparam idx Velocity vector index.
     */
    template<int idx>
    struct Moments{
        static constexpr int ma_Moments[Q]={1,
                                        3 * (CModulus<idx>)-4,
                                        Ci_x[idx],
                                        Ci_y[idx],
                                        Ci_x[idx]*Ci_x[idx]-Ci_y[idx]*Ci_y[idx]};
    };

    static constexpr std::array<int,Q*Q> MRTMatrix = GenerateMRTMatrix<Q,Moments>(); //!<MRT matrix array.

    /**
     * \brief Function that returns the diagonal entries of the MRT weight matrix for a given inverse relaxation time.
     * \param invtau Inverse of the relaxation time.
     */
    inline static std::vector<double> MRTWeights(const double& invtau){
        return {0,1,0,0,invtau};
    }

};

std::map<std::array<int8_t,D2Q5::D>,int> D2Q5::QMap = {{{0,0},0},
                                                      {{1,0},1},
                                                      {{-1,0},2},
                                                      {{0,1},3},
                                                      {{0,-1},4}};

/**
 * \brief 2D Stencil with diagonal directions.
 */
struct D2Q9:StencilBase {
    
    static constexpr int D = 2; //!<Number of cartesian directions.
    static constexpr int Q = 9; //!<Number of velocity directions.
    static constexpr double Cs2 = 1.0/3.0; //!<Speed of sound squared.
    
    static constexpr int Ci_x[Q] = {0, 1, -1, 0, 0, 1, -1, 1, -1}; //!<X component vector of velocity directions.
    static constexpr int Ci_y[Q] = {0, 0, 0, 1, -1, 1, -1, -1, 1}; //!<Y component vector of velocity directions.
    static constexpr int Ci_z[Q] = {0, 0, 0, 0, 0, 0, 0, 0, 0}; //!<Z component vector of velocity directions. 0 for 2D stencil.

    static std::map<std::array<int8_t,D>,int> QMap; //!<Map of cartesian vector to velocity index.

    enum{x = 0, y = 1, z = 2};
    
    /**
     * \brief Returns velocity direction vector depending on input d.
     * \param d Cartesian direction (x=0, y=1, z=2).
     * \return Vector of d component of velocity vectors.
     */
    inline static auto Ci_xyz(const int d) -> const int(&)[Q] {

        if (d == x) {
            return Ci_x;
        }
        else if (d == y) {
            return Ci_y;
        }
        return Ci_z;

    }
    static constexpr int Opposites[Q] = {0, 2, 1, 4, 3, 6, 5, 8, 7}; //!<Opposite vector at a given index.

    static constexpr double Weights[Q] = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0}; //!<Lattice weights.

    template<int idx>
    static constexpr int CModulus = Ci_x[idx] * Ci_x[idx] + Ci_y[idx] * Ci_y[idx]; //!<Computes the modulus of the velocity vector at a given index, used for the MRT weight calculation.

    /**
     * \brief Class containing an array of the MRT matrix moments for a given idx.
     * \tparam idx Velocity vector index.
     */
    template<int idx>
    struct Moments{
        static constexpr int ma_Moments[Q] = {1,
                                            3 * (CModulus<idx>) - 4,
                                            (9 * (CModulus<idx> * CModulus<idx>) - 21 * (CModulus<idx>) + 8) / 2,
                                            Ci_x[idx],
                                            (3 * (CModulus<idx>) - 5) * Ci_x[idx],
                                            Ci_y[idx],
                                            (3 * (CModulus<idx>) - 5) * Ci_y[idx],
                                            Ci_x[idx] * Ci_x[idx] - Ci_y[idx] * Ci_y[idx],
                                            Ci_x[idx] * Ci_y[idx]};
    };

    static constexpr std::array<int,Q*Q> MRTMatrix = GenerateMRTMatrix<Q,Moments>(); //!<MRT matrix.

    /**
     * \brief Function that returns the diagonal entries of the MRT weight matrix for a given inverse relaxation time.
     * \param invtau Inverse of the relaxation time.
     */
    inline static std::vector<double> MRTWeights(const double& invtau) {

        return {0, 1, 1, 0, 1, 0, 1, invtau, invtau};
        
    }
    
};

std::map<std::array<int8_t,D2Q9::D>,int> D2Q9::QMap = {{{0,0},0},
                                                      {{1,0},1},
                                                      {{-1,0},2},
                                                      {{0,1},3},
                                                      {{0,-1},4},
                                                      {{1,1},5},
                                                      {{-1,-1},6},
                                                      {{1,-1},7},
                                                      {{-1,1},8}};

/**
 * \brief Simple 3D stencil. Limited number of diagonal directions. Galilean invariance issues.
 */
struct D3Q15:StencilBase{

    static constexpr int D = 3; //!<Number of cartesian directions.
    static constexpr int Q = 15; //!<Number of velocity directions.
    static constexpr double Cs2 = 1.0/3.0; //!<Speed of sound squared.

    static constexpr int Ci_x[Q] = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1}; //!<X component vector of velocity directions.
    static constexpr int Ci_y[Q] = {0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1}; //!<Y component vector of velocity directions.
    static constexpr int Ci_z[Q] = {0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1}; //!<Z component vector of velocity directions.

    static std::map<std::array<int8_t,D>,int> QMap; //!<Map of cartesian vector to velocity index.

    /**
     * \brief Returns velocity direction vector depending on input d.
     * \param d Cartesian direction (x=0, y=1, z=2).
     * \return Vector of d component of velocity vectors.
     */
    inline static auto Ci_xyz(const int d) -> const int(&)[Q]{
        if (d==0) {
            return Ci_x;
        }
        else if (d==1) {
            return Ci_y;
        }
        return Ci_z;
    }

    static constexpr int Opposites[Q] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13}; //!<Opposite vector at a given index.
    static constexpr double Weights[Q] = {2.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0, 
				                          1.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,  1.0 / 72.0, 
				                          1.0 / 72.0,  1.0 / 72.0,  1.0 / 72.0,  1.0 / 72.0, 
				                          1.0 / 72.0,  1.0 / 72.0,  1.0 / 72.0}; //!<Lattice weights.

    template<int idx>
    static const int CModulus = Ci_x[idx] * Ci_x[idx] + Ci_y[idx] * Ci_y[idx] + Ci_z[idx] * Ci_z[idx]; //!<Computes the modulus of the velocity vector at a given index, used for the MRT weight calculation.

    /**
     * \brief Class containing an array of the MRT matrix moments for a given idx.
     * \tparam idx Velocity vector index.
     */
    template<int idx>
    struct Moments{
        static constexpr int ma_Moments[Q] = {1,
                                        (CModulus<idx>) - 2,
                                        (15 * (CModulus<idx> * CModulus<idx>) - 55 * (CModulus<idx>) + 32) / 2,
                                        Ci_x[idx],
                                        (5 * (CModulus<idx>) - 13) * Ci_x[idx],
                                        Ci_y[idx],
                                        (5 * (CModulus<idx>) - 13) * Ci_y[idx],
                                        Ci_z[idx],
                                        (5 * (CModulus<idx>) - 13) * Ci_z[idx],
                                        3 * Ci_x[idx] * Ci_x[idx] - CModulus<idx>,
                                        Ci_y[idx] * Ci_y[idx] - Ci_z[idx] * Ci_z[idx],
                                        Ci_x[idx] * Ci_y[idx],
                                        Ci_y[idx] * Ci_z[idx],
                                        Ci_x[idx] * Ci_z[idx],
                                        Ci_x[idx] * Ci_y[idx] * Ci_z[idx]};
    };

    static constexpr std::array<int,Q*Q> MRTMatrix = GenerateMRTMatrix<Q,Moments>(); //!<MRT matrix.

    /**
     * \brief Function that returns the diagonal entries of the MRT weight matrix for a given inverse relaxation time.
     * \param invtau Inverse of the relaxation time.
     */
    inline static std::vector<double> MRTWeights(const double& invtau) {

        return {0, 1, 1, 0, 1, 0, 1, 0, 1, invtau, invtau, invtau, invtau, invtau, 1, 1};

    }
};

std::map<std::array<int8_t,D3Q15::D>,int> D3Q15::QMap = {{{0,0,0},0},
                                                       {{1,0,0},1},
                                                       {{-1,0,0},2},
                                                       {{0,1,0},3},
                                                       {{0,-1,0},4},
                                                       {{0,0,1},5},
                                                       {{0,0,-1},6},
                                                       {{1,1,1},7},
                                                       {{-1,-1,-1},8},
                                                       {{1,1,-1},9},
                                                       {{-1,-1,1},10},
                                                       {{1,-1,1},11},
                                                       {{-1,1,-1},12},
                                                       {{-1,1,1},13},
                                                       {{1,-1,-1},14}};

/**
 * \brief Most commonly used 3D stencil.
 */
struct D3Q19:StencilBase{

    static constexpr int D = 3; //!<Number of cartesian directions.
    static constexpr int Q = 19; //!<Number of velocity directions.
    static constexpr double Cs2 = 1.0/3.0; //!<Speed of sound squared.

    static constexpr int Ci_x[Q] = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1}; //!<X component vector of velocity directions.
    static constexpr int Ci_y[Q] = {0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 1, -1, 1, -1, 0, 0, 0, 0}; //!<Y component vector of velocity directions.
    static constexpr int Ci_z[Q] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1}; //!<Z component vector of velocity directions.

    static std::map<std::array<int8_t,D>,int> QMap; //!<Map of cartesian vector to velocity index.

    /**
     * \brief Returns velocity direction vector depending on input d.
     * \param d Cartesian direction (x=0, y=1, z=2).
     * \return Vector of d component of velocity vectors.
     */
    inline static auto Ci_xyz(const int d) -> const int(&)[Q]{
        if (d==0) {
            return Ci_x;
        }
        else if (d==1) {
            return Ci_y;
        }
        return Ci_z;
    }

    static constexpr int Opposites[Q] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17}; //!<Opposite vector at a given index.
    static constexpr double Weights[Q] = {1.0 / 3.0,  1.0 / 18.0,  1.0 / 18.0,  1.0 / 18.0, 
				                          1.0 / 18.0,  1.0 / 18.0,  1.0 / 18.0,  1.0 / 36.0, 
				                          1.0 / 36.0,  1.0 / 36.0,  1.0 / 36.0,  1.0 / 36.0, 
				                          1.0 / 36.0,  1.0 / 36.0,  1.0 / 36.0,  1.0 / 36.0, 
				                          1.0 / 36.0,  1.0 / 36.0,  1.0 / 36.0}; //!<Lattie weights.


    template<int idx>
    static const int CModulus = Ci_x[idx] * Ci_x[idx] + Ci_y[idx] * Ci_y[idx] + Ci_z[idx] * Ci_z[idx]; //!<Computes the modulus of the velocity vector at a given index, used for the MRT weight calculation.

    /**
     * \brief Class containing an array of the MRT matrix moments for a given idx.
     * \tparam idx Velocity vector index.
     */
    template<int idx>
    struct Moments{
        static constexpr int ma_Moments[Q] = {1,
                                        19 * (CModulus<idx>) - 30,
                                        (21 * (CModulus<idx> * CModulus<idx>) - 53 * (CModulus<idx>) + 24) / 2,
                                        Ci_x[idx],
                                        (5 * (CModulus<idx>) - 9) * Ci_x[idx],
                                        Ci_y[idx],
                                        (5 * (CModulus<idx>) - 9) * Ci_y[idx],
                                        Ci_z[idx],
                                        (5 * (CModulus<idx>) - 9) * Ci_z[idx],
                                        3 * Ci_x[idx] * Ci_x[idx] - CModulus<idx>,
                                        (3 * CModulus<idx> - 5) * (3 * Ci_x[idx] * Ci_x[idx] - CModulus<idx>),
                                        Ci_y[idx] * Ci_y[idx] - Ci_z[idx] * Ci_z[idx],
                                        (3 * CModulus<idx> - 5) * (Ci_y[idx] * Ci_y[idx] - Ci_z[idx] * Ci_z[idx]),
                                        Ci_x[idx] * Ci_y[idx],
                                        Ci_y[idx] * Ci_z[idx],
                                        Ci_x[idx] * Ci_z[idx],
                                        Ci_x[idx] * (Ci_y[idx] * Ci_y[idx] - Ci_z[idx] * Ci_z[idx]),
                                        Ci_y[idx] * (Ci_z[idx] * Ci_z[idx] - Ci_x[idx] * Ci_x[idx]),
                                        Ci_z[idx] * (Ci_x[idx] * Ci_x[idx] - Ci_y[idx] * Ci_y[idx])};
    };

    static constexpr std::array<int,Q*Q> MRTMatrix = GenerateMRTMatrix<Q,Moments>(); //!<MRT matrix.

    /**
     * \brief Function that returns the diagonal entries of the MRT weight matrix for a given inverse relaxation time.
     * \param invtau Inverse of the relaxation time.
     */
    inline static std::vector<double> MRTWeights(const double& invtau) {

        return {0, 1, 1, 0, 1, 0, 1, 0, 1, invtau, 1, invtau, 1, invtau, invtau, invtau, 1, 1, 1};

    }
};

std::map<std::array<int8_t,D3Q19::D>,int> D3Q19::QMap = {{{0,0,0},0},
                                                       {{1,0,0},1},
                                                       {{-1,0,0},2},
                                                       {{0,1,0},3},
                                                       {{0,-1,0},4},
                                                       {{0,0,1},5},
                                                       {{0,0,-1},6},
                                                       {{1,1,0},7},
                                                       {{-1,-1,0},8},
                                                       {{1,-1,0},9},
                                                       {{-1,1,0},10},
                                                       {{0,1,1},11},
                                                       {{0,-1,-1},12},
                                                       {{0,1,-1},13},
                                                       {{0,-1,1},14},
                                                       {{1,0,1},15},
                                                       {{-1,0,-1},16},
                                                       {{1,0,-1},17},
                                                       {{-1,0,1},18}};

/**
 * Higher isotropy 3D stencil.
 */
struct D3Q27:StencilBase{

    static constexpr int D = 3; //!<Number of cartesian directions.
    static constexpr int Q = 27; //!<Number of velocity directions.
    static constexpr double Cs2 = 1.0/3.0; //!<Speed of sound squared.

    static constexpr int Ci_x[Q] = {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1}; //!<X component vector of velocity directions.
    static constexpr int Ci_y[Q] = {0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1}; //!<Y component vector of velocity directions.
    static constexpr int Ci_z[Q] = {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1}; //!<Z component vector of velocity directions.

    static std::map<std::array<int8_t,D>,int> QMap; //!<Map of cartesian vector to velocity index.

    /**
     * \brief Returns velocity direction vector depending on input d.
     * \param d Cartesian direction (x=0, y=1, z=2).
     * \return Vector of d component of velocity vectors.
     */
    inline static auto Ci_xyz(const int d) -> const int(&)[Q]{
        if (d==0) {
            return Ci_x;
        }
        else if (d==1) {
            return Ci_y;
        }
        return Ci_z;
    }

    static constexpr int Opposites[Q] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25}; //!<Opposite vector at a given index.
    static constexpr double Weights[Q] = {8.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0, 
				                          2.0 / 27.0,  2.0 / 27.0,  2.0 / 27.0,  1.0 / 54.0, 
				                          1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0, 
				                          1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0, 
				                          1.0 / 54.0,  1.0 / 54.0,  1.0 / 54.0, }; //!<Lattice weights.

    template<int idx>
    static const int CModulus = Ci_x[idx] * Ci_x[idx] + Ci_y[idx] * Ci_y[idx] + Ci_z[idx] * Ci_z[idx]; //!<Computes the modulus of the velocity vector at a given index, used for the MRT weight calculation.

};

std::map<std::array<int8_t,D3Q27::D>,int> D3Q27::QMap = {{{0,0,0},0},
                                                        {{1,0,0},1},
                                                        {{-1,0,0},2},
                                                        {{0,1,0},3},
                                                        {{0,-1,0},4},
                                                        {{0,0,1},5},
                                                        {{0,0,-1},6},
                                                        {{1,1,0},7},
                                                        {{-1,-1,0},8},
                                                        {{1,-1,0},9},
                                                        {{-1,1,0},10},
                                                        {{0,1,1},11},
                                                        {{0,-1,-1},12},
                                                        {{0,1,-1},13},
                                                        {{0,-1,1},14},
                                                        {{1,0,1},15},
                                                        {{-1,0,-1},16},
                                                        {{1,0,-1},17},
                                                        {{-1,0,1},18},
                                                        {{1,1,1},19},
                                                        {{-1,-1,-1},20},
                                                        {{1,1,-1},21},
                                                        {{-1,-1,1},22},
                                                        {{1,-1,1},23},
                                                        {{-1,1,-1},24},
                                                        {{-1,1,1},25},
                                                        {{1,-1,-1},26}};
