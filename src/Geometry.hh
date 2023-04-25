#ifndef GEOMETRY_HEADER
#define GEOMETRY_HEADER
#include<map>
#include "Service.hh"
 
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
 */
class Geometry{
    public:

        /**
         * \brief Initialise the geometry vector for the system based on some conditions in this function.
         */
        void initialise();

        /**
         * \brief Returns true if the current lattice point lies on a periodic boundary.
         * \param k Index of current lattice point.
         */
        bool isPeriodic(int k);

        /**
         * \brief Returns true if the current lattice point lies on a solid boundary.
         * \param k Index of current lattice point.
         */
        bool isSolid(int k);

    private:

        enum{Periodic=3,Solid=2,Wall=1,Fluid=0}; //!< IDs for each boundary type.

        std::vector<int> geometry; //!< Vector containing geometry information.


};

/**
 * \details This function evaluates an (unnecessarily complex for now) if statement to determine if we are on the
 *          edge of the simulation domain. Returns true if current lattice point is a periodic boundary.
 */
bool Geometry::isPeriodic(int k){

    if (((k%(LZ-1)==0||
         (k-1)%(LZ-1)==0)&&LZ>1)||
        ((k/(LZ)%(LY-1)==0||
         ((k)/(LZ)-1)%(LY-1)==0)&&LY>1)||
        (((k/(LZ)/(LY))%(LXdiv-1)==0||
          ((k)/(LZ)/(LY)-1)%(LXdiv-1)==0||
          ((k)/(LZ)/(LY)-1)-1<0||
          ((k)/(LZ)/(LY)-1)+1>LXdiv-1)&&
         LXdiv>1)||
        LZ==1||
        LY==1||
        LXdiv==1) return true; //These conditions just check whether the lattice point is on the edge of the 
                            //simulation domain
    
    else return false;
}

/**
 * \details This function evaluates an if statement to determine if we are on the solid boundary. This boundary
 *          is currently determined by this function, but in the future it will be chosen elsewhere so that no
 *          modification to the source file is needed. Returns true if current lattice point is a solid.
 */
bool Geometry::isSolid(int k){

    int yAtCurrentk=computeY(k);

    if (yAtCurrentk<=1||yAtCurrentk>=LY-2) return true; //Change this condition to control where the solid is
    
    else return false;

}

#endif