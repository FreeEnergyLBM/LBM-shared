#pragma once
#include<map>
#include "Parameters.hh"
#include "Service.hh"
#include "Lattice.hh"
#include "Global.hh"

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
         * \brief Returns true if the current TLattice point lies on a periodic boundary.
         * \param k Index of current TLattice point.
         * \return True if TLattice point lies on a periodic boundary.
         */
        static inline bool isPeriodic(int k);

        /**
         * \brief Returns true if the current TLattice point lies on a solid boundary.
         * \param k Index of current TLattice point.
         * \return True if TLattice point lies on a solid
         */
        static inline bool isSolid(int k);

    private:

        enum{Periodic = 3, Solid = 2, Wall = 1, Fluid = 0}; //!<IDs for each boundary type.

};

/**
 * \details This function evaluates an (unnecessarily complex for now) if statement to determine if we are on the
 *          edge of the simulation domain. Returns true if current TLattice point is a periodic boundary.
 */
template<class TLattice>
inline bool Geometry<TLattice>::isPeriodic(int k) {
    
    int yAtCurrentk = computeY(TLattice::LY, TLattice::LZ, k);
    int zAtCurrentk = computeZ(TLattice::LY, TLattice::LZ, k);
    int xAtCurrentk = computeX(TLattice::LY, TLattice::LZ, k);

    if(TLattice::LZ <= 1 || TLattice::LY <= 1 || TLattice::LXdiv <= 1) return true; //If simulation is 2D
    else if (zAtCurrentk == 0 ||
         zAtCurrentk == TLattice::LZ-1) return true; //Edges in Z direction

    else if (yAtCurrentk == 0 ||
         yAtCurrentk == TLattice::LY-1) return true; //Edges in Y direction
        
    else if (xAtCurrentk == 0 ||
         xAtCurrentk == TLattice::LXdiv-1) return true; //Edges in X direction
    
    return false;
    
    /*
    if (((TLattice::LZ>1&&(k%(TLattice::LZ-1)==0||
         (k-1)%(TLattice::LZ-1)==0)))||
        (TLattice::LY>1&&(k/(TLattice::LZ)%(TLattice::LY-1)==0||
         ((k)/(TLattice::LZ)-1)%(TLattice::LY-1)==0))||
        (TLattice::LXdiv>1&&((k/(TLattice::LZ)/(TLattice::LY))%(TLattice::LXdiv-1)==0||
          ((k)/(TLattice::LZ)/(TLattice::LY)-1)%(TLattice::LXdiv-1)==0||
          ((k)/(TLattice::LZ)/(TLattice::LY)-1)-1<0||
          ((k)/(TLattice::LZ)/(TLattice::LY)-1)+1>TLattice::LXdiv-1))||
        TLattice::LZ==1||
        TLattice::LY==1||
        TLattice::LXdiv==1) return true; //These conditions just check whether the TLattice point is on the edge of the 
                            //simulation domain
    
    else return false;
    */
}

/**
 * \details This function evaluates an if statement to determine if we are on the solid boundary. This boundary
 *          is currently determined by this function, but in the future it will be chosen elsewhere so that no
 *          modification to the source file is needed. Returns true if current TLattice point is a solid.
 */
template<class TLattice>
inline bool Geometry<TLattice>::isSolid(int k) {

    return SolidLabels<>::get<TLattice>(k);
    //int yAtCurrentk = computeY(TLattice::LY, TLattice::LZ, k);
    //int zAtCurrentk = computeZ(TLattice::LY, TLattice::LZ, k);
    //int xAtCurrentk = computeX(TLattice::LY, TLattice::LZ, k);
    
    //if (yAtCurrentk <= 1 || yAtCurrentk >= TLattice::LY - 2 || xAtCurrentk <= 1 || xAtCurrentk >= TLattice::LXdiv - 2 || zAtCurrentk <= 1 || zAtCurrentk >= TLattice::LZ - 2 ) return true; //Change this condition to control where the solid is
    
    //if (yAtCurrentk <= 1 || yAtCurrentk >= TLattice::LY - 2 ) return true; //Change this condition to control where the solid is
    //else return false;

}