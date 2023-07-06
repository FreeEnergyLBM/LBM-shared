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
 * This class contains functions to initialise the chosen geometry and determine whether the current lattice
 * point is a solid or a periodic boundary.
 * \tparam Empty template parameter, used to ensure GETPROPERTIES is defined before it is needed in this function.
 */

template<class lattice>
class Geometry {
    public:

        /**
         * \brief Returns true if the current lattice point lies on a periodic boundary.
         * \param k Index of current lattice point.
         * \return True if lattice point lies on a periodic boundary.
         */
        static inline bool isPeriodic(int k);

        /**
         * \brief Returns true if the current lattice point lies on a solid boundary.
         * \param k Index of current lattice point.
         * \return True if lattice point lies on a solid
         */
        static inline bool isSolid(int k);

    private:

        enum{Periodic = 3, Solid = 2, Wall = 1, Fluid = 0}; //!<IDs for each boundary type.

};

/**
 * \details This function evaluates an (unnecessarily complex for now) if statement to determine if we are on the
 *          edge of the simulation domain. Returns true if current lattice point is a periodic boundary.
 */
template<typename lattice>
inline bool Geometry<lattice>::isPeriodic(int k) {
    
    int yAtCurrentk = computeY(lattice::LY, lattice::LZ, k);
    int zAtCurrentk = computeZ(lattice::LY, lattice::LZ, k);
    int xAtCurrentk = computeX(lattice::LY, lattice::LZ, k);

    if(lattice::LZ <= 1 || lattice::LY <= 1 || lattice::LXdiv <= 1) return true; //If simulation is 2D
    else if (zAtCurrentk == 0 ||
         zAtCurrentk == lattice::LZ-1) return true; //Edges in Z direction

    else if (yAtCurrentk == 0 ||
         yAtCurrentk == lattice::LY-1) return true; //Edges in Y direction
        
    else if (xAtCurrentk == 0 ||
         xAtCurrentk == lattice::LXdiv-1) return true; //Edges in X direction
    
    return false;
    
    /*
    if (((lattice::LZ>1&&(k%(lattice::LZ-1)==0||
         (k-1)%(lattice::LZ-1)==0)))||
        (lattice::LY>1&&(k/(lattice::LZ)%(lattice::LY-1)==0||
         ((k)/(lattice::LZ)-1)%(lattice::LY-1)==0))||
        (lattice::LXdiv>1&&((k/(lattice::LZ)/(lattice::LY))%(lattice::LXdiv-1)==0||
          ((k)/(lattice::LZ)/(lattice::LY)-1)%(lattice::LXdiv-1)==0||
          ((k)/(lattice::LZ)/(lattice::LY)-1)-1<0||
          ((k)/(lattice::LZ)/(lattice::LY)-1)+1>lattice::LXdiv-1))||
        lattice::LZ==1||
        lattice::LY==1||
        lattice::LXdiv==1) return true; //These conditions just check whether the lattice point is on the edge of the 
                            //simulation domain
    
    else return false;
    */
}

/**
 * \details This function evaluates an if statement to determine if we are on the solid boundary. This boundary
 *          is currently determined by this function, but in the future it will be chosen elsewhere so that no
 *          modification to the source file is needed. Returns true if current lattice point is a solid.
 */
template<typename lattice>
inline bool Geometry<lattice>::isSolid(int k) {

    return SolidLabels<>::get<lattice>(k);
    //int yAtCurrentk = computeY(lattice::LY, lattice::LZ, k);
    //int zAtCurrentk = computeZ(lattice::LY, lattice::LZ, k);
    //int xAtCurrentk = computeX(lattice::LY, lattice::LZ, k);
    
    //if (yAtCurrentk <= 1 || yAtCurrentk >= lattice::LY - 2 || xAtCurrentk <= 1 || xAtCurrentk >= lattice::LXdiv - 2 || zAtCurrentk <= 1 || zAtCurrentk >= lattice::LZ - 2 ) return true; //Change this condition to control where the solid is
    
    //if (yAtCurrentk <= 1 || yAtCurrentk >= lattice::LY - 2 ) return true; //Change this condition to control where the solid is
    //else return false;

}