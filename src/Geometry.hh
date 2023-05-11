#pragma once
#include<map>
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
template<typename placeholder = void>
class GeometryTemplate {
    public:

        /**
         * \brief Initialise the geometry vector for the system based on some conditions in this function.
         */
        inline void initialise();

        /**
         * \brief Returns true if the current lattice point lies on a periodic boundary.
         * \param k Index of current lattice point.
         * \return True if lattice point lies on a periodic boundary.
         */
        inline bool isPeriodic(int k);

        /**
         * \brief Returns true if the current lattice point lies on a solid boundary.
         * \param k Index of current lattice point.
         * \return True if lattice point lies on a solid
         */
        inline bool isSolid(int k);

    private:

        enum{Periodic = 3, Solid = 2, Wall = 1, Fluid = 0}; //!<IDs for each boundary type.

        std::vector<int> geometry; //!<Vector containing geometry information.


};

/**
 * \details This function evaluates an (unnecessarily complex for now) if statement to determine if we are on the
 *          edge of the simulation domain. Returns true if current lattice point is a periodic boundary.
 */
template<typename placeholder>
inline bool GeometryTemplate<placeholder>::isPeriodic(int k) {

    if(GETPROPERTIES().m_LZ <= 1 || GETPROPERTIES().m_LY <= 1 || GETPROPERTIES().m_LXdiv <= 1) return true; //If simulation is 2D
    else if (k % (GETPROPERTIES().m_LZ - 1) == 0 ||
         (k - 1) % (GETPROPERTIES().m_LZ - 1) == 0) return true; //Edges in Z direction

    else if (k / (GETPROPERTIES().m_LZ) % (GETPROPERTIES().m_LY - 1) == 0 ||
         ((k) / (GETPROPERTIES().m_LZ) - 1) % (GETPROPERTIES().m_LY - 1) == 0) return true; //Edges in Y direction
        
    else if ((k / (GETPROPERTIES().m_LZ) / (GETPROPERTIES().m_LY)) % (GETPROPERTIES().m_LXdiv - 1) == 0 ||
          ((k) / (GETPROPERTIES().m_LZ) / (GETPROPERTIES().m_LY) - 1) % (GETPROPERTIES().m_LXdiv - 1) == 0 ||
          ((k) / (GETPROPERTIES().m_LZ) / (GETPROPERTIES().m_LY) - 1) - 1 < 0 ||
          ((k) / (GETPROPERTIES().m_LZ) / (GETPROPERTIES().m_LY) - 1) + 1 > GETPROPERTIES().m_LXdiv - 1) return true; //Edges in X direction
    
    return false;

}

/**
 * \details This function evaluates an if statement to determine if we are on the solid boundary. This boundary
 *          is currently determined by this function, but in the future it will be chosen elsewhere so that no
 *          modification to the source file is needed. Returns true if current lattice point is a solid.
 */
template<typename placeholder>
inline bool GeometryTemplate<placeholder>::isSolid(int k) {

    int yAtCurrentk = computeY(GETPROPERTIES().m_LY, GETPROPERTIES().m_LZ, k);

    if (yAtCurrentk <= 1 || yAtCurrentk >= GETPROPERTIES().m_LY - 2) return true; //Change this condition to control where the solid is
    
    else return false;

}

typedef GeometryTemplate<> Geometry; //!<Typedef to avoid needing empty angle brackets