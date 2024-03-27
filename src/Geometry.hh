#pragma once
#include<map>
#include "Parameters.hh"
#include "Service.hh"
#include "Data.hh"

/**
 * \file Geometry.hh
 * \brief This is not fully implemented but will contain information relevant to boundary conditions
 * This contains a geometry class that will contain information about the solid geometry and boundary locations
 * in the simulation.
 */

/**
 * \brief Geometry contains functions to initialise and access the chosen geometry.
 * This class contains functions to initialise the chosen geometry and determine whether the current TLattice
 * point is a solid or a periodic boundary.
 * \tparam Empty template parameter, used to ensure GETPROPERTIES is defined before it is needed in this function.
 */
template<class TLattice>
class Geometry {
    public:

        /**
         * \brief Returns true if the current TLattice point lies on a solid boundary.
         * \param k Index of current TLattice point.
         * \return True if TLattice point lies on a solid
         */
        static inline bool isBoundary(int k);

        template<class TStencil>
        static inline std::array<int8_t,TLattice::NDIM> findNormal(int (*condition)(const int),const std::vector<int>& neighbors,const std::vector<int> fluidvals,int k);

        static inline bool isCorner(int (*condition)(const int),const std::array<int8_t,TLattice::NDIM>& normal,const std::vector<int> fluidvals,int k);

        static inline bool isCorner(int k);

        static inline bool isBulkSolid(int k,int (*condition)(const int),const std::vector<int>& neighbors,std::vector<int> fluidvals);

        static inline bool isBulkSolid(int k);

        static inline void initialiseBoundaries(int (*condition)(const int),int fluidval = 0);

        static inline void initialiseBoundaries(int (*condition)(const int),std::vector<int> fluidvals);

        static inline int& getBoundaryType(int k);

        static std::vector<int> mFluidVals;

    private:

        enum
        {
            BulkSolid = -1,
            Fluid = 0,
            Wall = 1,
            Wall2 = 2,
            InletWall = 3,
            OutletWall = 4,
            HumidityInterface = 5,
            RefillNode = 6
        };

        static inline bool apply(int k, int (*condition)(const int),std::vector<int> fluidvals) {
            
            for (int i : fluidvals){
                if(condition(k) == i) return false;
            }
            return true;
        }

};

template<class TLattice>
std::vector<int> Geometry<TLattice>::mFluidVals;

template<class TLattice>
inline bool Geometry<TLattice>::isBulkSolid(int k,int (*condition)(const int),const std::vector<int>& neighbors,std::vector<int> fluidvals) {

    mFluidVals=fluidvals;

    using Stencil = std::conditional_t<TLattice::NDIM == 1, D1Q3, std::conditional_t<TLattice::NDIM == 2, D2Q9, D3Q27>>;

    bool bulksolid = true;
        
    for (int idx = 0; idx < Stencil::Q; idx++) {
        for (int i : fluidvals){
            if(condition(neighbors[k*Stencil::Q+idx]) == i) bulksolid=false;
        }
    }

    return bulksolid;

}

template<class TLattice>
inline void Geometry<TLattice>::initialiseBoundaries(int (*condition)(const int),const std::vector<int> fluidvals) {

    using Stencil = std::conditional_t<TLattice::NDIM == 1, D1Q3, std::conditional_t<TLattice::NDIM == 2, D2Q9, D3Q27>>;

    using data = Data_Base<TLattice, Stencil>;

    std::vector<int>& neighbors = data::getInstance().getNeighbors();

    for(int k = TLattice::HaloSize; k < TLattice::N-TLattice::HaloSize; k++) {

        int solidval;
        
        if (isBulkSolid(k,condition,neighbors,fluidvals)) {
            solidval = -1;
        }
        else {
            solidval = condition(k);
        }

        std::array<int8_t,TLattice::NDIM> normal=findNormal<Stencil>(condition,neighbors,fluidvals,k);

        Boundary<TLattice::NDIM> boundaryk = {solidval,isCorner(condition,normal,fluidvals,k),normal};
      
        BoundaryLabels<TLattice::NDIM>::template initialise<TLattice>(boundaryk,k);

    }

    neighbors.resize(0);

}

template<class TLattice>
inline void Geometry<TLattice>::initialiseBoundaries(int (*condition)(const int), int fluidval) {

    initialiseBoundaries(condition,std::vector<int>{fluidval});

}

/**
 * \details This function evaluates an if statement to determine if we are on the solid boundary. This boundary
 *          is currently determined by this function, but in the future it will be chosen elsewhere so that no
 *          modification to the source file is needed. Returns true if current TLattice point is a solid.
 */
template<class TLattice>
inline bool Geometry<TLattice>::isBoundary(int k) {
    for (int i : mFluidVals){
        if(BoundaryLabels<TLattice::NDIM>::template get<TLattice>(k).Id == i) return false;
    }
    return true;

}

template<class TLattice>
inline bool Geometry<TLattice>::isCorner(int k) {

    return BoundaryLabels<TLattice::NDIM>::template get<TLattice>(k).IsCorner;

}

template<class TLattice>
inline bool Geometry<TLattice>::isBulkSolid(int k) {

    return (BoundaryLabels<TLattice::NDIM>::template get<TLattice>(k).Id==-1);

}

template<class TLattice>
template<class TStencil>
inline std::array<int8_t,TLattice::NDIM> Geometry<TLattice>::findNormal(int (*condition)(const int), const std::vector<int>& neighbors,const std::vector<int> fluidvals,int k) {

    std::array<int8_t,TLattice::NDIM> normal = {};
    
    //if (!apply(k,condition,fluidvals)) return normal;

    std::vector<int> sum(TLattice::NDIM,0);

    for (int idx = 0; idx < TStencil::Q; idx++){

        if(!apply(neighbors[k*TStencil::Q+idx],condition,fluidvals)&&condition(k)!=condition(neighbors[k*TStencil::Q+idx])) {
            sum[0]+=TStencil::Ci_xyz(0)[idx];
            if constexpr (TLattice::NDIM>1) sum[1]+=TStencil::Ci_xyz(1)[idx];
            if constexpr (TLattice::NDIM>2) sum[2]+=TStencil::Ci_xyz(2)[idx];
        }

    }

    if constexpr (TLattice::NDIM>1){
        if(abs(sum[0])>abs(sum[1])){
            if constexpr (TLattice::NDIM>2){
                if (abs(sum[0])>abs(sum[2])){
                    sum[2]=0;
                }
                else if (abs(sum[0])<abs(sum[2])){
                    sum[0]=0;
                }
                if (abs(sum[1])>abs(sum[2])){
                    sum[2]=0;
                }
                else if (abs(sum[1])<abs(sum[2])){
                    sum[1]=0;
                }
            }
            sum[1]=0;
        }
        else if (abs(sum[0])<abs(sum[1])) {
            if constexpr (TLattice::NDIM>2){
                if (abs(sum[0])>abs(sum[2])){
                    sum[2]=0;
                }
                else if (abs(sum[0])<abs(sum[2])){
                    sum[0]=0;
                }
                if (abs(sum[1])>abs(sum[2])){
                    sum[2]=0;
                }
                else if (abs(sum[1])<abs(sum[2])){
                    sum[1]=0;
                }
            }
            sum[0]=0;
        }
    }

    for (int xyz = 0; xyz < TLattice::NDIM; xyz++){

        normal[xyz] = (int8_t) ((sum[xyz]==0) ? 0 : ((sum[xyz] > 0) ? 1 : -1));

    }

    return normal;

}

template<class TLattice>
inline bool Geometry<TLattice>::isCorner(int (*condition)(const int),const std::array<int8_t,TLattice::NDIM>& normal,const std::vector<int> fluidvals,int k) { // DOesnt

    if (!apply(k,condition,fluidvals)) return false;

    int normalsum = 0;

    for (int i=0; i < (int)normal.size(); i++) { 
        normalsum += abs((int)normal[i]);
    }

    return normalsum > 1;

}

template<class TLattice>
inline int& Geometry<TLattice>::getBoundaryType(int k) {
    //static auto& test = BoundaryLabels<TLattice::NDIM>::template get<TLattice>();
    return BoundaryLabels<TLattice::NDIM>::template get<TLattice>(k).Id;

}

