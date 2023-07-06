#pragma once
#include "Lattice.hh"
#include "Geometry.hh"
#include "Parameters.hh"
#include "Service.hh"

/**
 * \file Data.hh
 * \brief Contains data class that will control how how data is accessed and how stremaing happens.
 * These classes will contain derived distribution classes that will determine memory allocation for the
 * distribution arrays and the streaming indices. The Data class also contains a vector of neighbors for each 
 * lattice point in every direction but it might be faster to just recalculate every time this is needed.
 * "Neighbor" refers to the lattice point adjacent to any given latice point in a chosen discrete direction
 * Periodic boundaries work by setting the neighbor of each lattice point at the edges of the domain so that the
 * top of the domain connects to the bottom, the left connects to the right etc.
 */


/**
 * \brief The Data_Base class provides information for the layout of neighboring lattice points and communication
          for non-distribution parameters.
 * This class takes a stencil as a template argument, as the velocity discretisation information and weights is 
 * needed. The class has public functions for communication, generating neighbors based on the
 * stencil and a function to get a reference to the vector of neighbor indices.
 * \tparam stencil Velocity Stencil of class using this data type.
 */
template<class lattice, class stencil>
class Data_Base{
    public:
        static Data_Base<lattice,stencil>& getInstance() {
            static Data_Base<lattice,stencil> instance;
            return instance;
        }
    private:
        //ADD VIRTUAL TO THIS
        /**
         * \brief This function returns the neighbor at the current lattice point k in the direction Q.
         * \param k Index of current lattice point.
         * \param Q Discrete velocity direction (e.g. 0-8 for D2Q9).
         * \return Lattice index of neighboring lattice point in chosen direction.
         */
        inline int getOneNeighbor(const int k, const int Q);

        /**
         * \brief This function returns the neighbor at the current lattice point k in the direction Q given that
         *        the lattice point lies on a periodic boundary.
         * \param k Index of current lattice point.
         * \param Q Discrete velocity direction (e.g. 0-8 for D2Q9).
         * \return Lattice index of neighboring lattice point in chosen direction given that the current point is on a periodic boundary.
         */
        inline int getOneNeighborPeriodic(const int k, const int Q);

        std::array<int,stencil::Q> OppositeOffset; //!<Opposite lattice point offset in each direction.
        
        Geometry<lattice> m_Geometry; //!<Class containing geometry information (for periodic boundaries).

        std::vector<int> mv_Neighbors; //!<Vector containing neighbor information.

        enum{ x = 0, y = 1, z = 2 }; //!<Indices corresponding to x, y, z.
        
        template<class, class>
        friend class Data1; //Data1 can access private members of the base class (will need to add new data
                            //types here but I will probably change how this works).

        /**
         * \brief The constructor for the class.
         * This constructor will calculate the opposite points at each index Q,
         * allocate memory for neighbors and fill the array of neighbors.
         */
        Data_Base() //Construct distribution
        {

            for(int idx = 0; idx <stencil::Q; idx++) {

                OppositeOffset[idx] = stencil::Ci_xyz(x)[idx] * lattice::LZ * lattice::LY + stencil::Ci_xyz(y)[idx] * lattice::LZ + stencil::Ci_xyz(z)[idx];

            }

            mv_Neighbors.resize(stencil::Q * lattice::N); //Allocate memory for neighbors array
            
            generateNeighbors(); //Fill neighbors array
            
        }
        
    public:

        Data_Base(Data_Base<lattice, stencil> const& other) = delete;
        void operator=(Data_Base<lattice, stencil> const& other) = delete;

        using Stencil = stencil; //!<Typedef so type info can be accessed from outside the class.

        /**
         * \brief This function communicates a chosen parameter (halo regions are exchanged with neighboring
         *        processors).
         * \param obj Object of chosen parameter.
         * \tparam parameter type of object to be communicated.
         */
        template<class parameter>
        inline void communicate(parameter& obj);

        /**
         * \brief Function to fill neighbor array with neighbor information.
         */
        inline void generateNeighbors();

        /**
         * \brief Returns the neighbor array
         */
        inline std::vector<int>& getNeighbors();
        inline int& getNeighbor(int k,int q){return mv_Neighbors[k*stencil::Q+q];};
        const inline std::vector<int>& getNeighbors() const;

};

/**
 * \details The function calls the communicate(obj) function from the parallel template class chosen when the
 *          class is used. This will perform the necessary communications so gradients etc. can be calculated
 *          across parallel regions.
 */
template<class lattice, class stencil>
template<class parameter>
inline void Data_Base<lattice,stencil>::communicate(parameter& obj) { //Not used in this data type

    //static_assert(is_base_of_template<Parameter,parameter>::value,"ERROR: The object passed to this function cannot be communicated.");
    lattice::communicate(obj);

}

/**
 * \details The function returns a reference to a vector containing the neighboring lattice point for every point
 *          in every direction in the stencil.
 */
template<class lattice, class stencil>
inline std::vector<int>& Data_Base<lattice, stencil>::getNeighbors() {
    
    return mv_Neighbors;

}

template<class lattice, class stencil>
const inline std::vector<int>& Data_Base<lattice, stencil>::getNeighbors() const {
    
    return mv_Neighbors;

}

/**
 * \details The neighbor of the current lattice point is calculated from the current lattice point + the offset
 *          in each direction, which is precomputed in the constructor and stored in a vector.
 */
template<class lattice, class stencil>
inline int Data_Base<lattice,stencil>::getOneNeighbor(const int k, const int Q) {
    
    return k + OppositeOffset[Q]; //The neighbor is the lattice point plus the opposite offset in direction Q
        
}

/**
 * \details If we are at the first lattice point, some of the adjacent points will be on the complete opposite
 *          side of the lattice so we must account for this. The function will work out if we are on the edge of
 *          of the simulation in the x, y and z directions and apply offsets in each case.
 */
template<class lattice, class stencil>
inline int Data_Base<lattice,stencil>::getOneNeighborPeriodic(const int k, const int Q) {

    int neighbor = 0;

    if(lattice::LZ> 1) {
        int localz=computeZ(lattice::LY,lattice::LZ,k);
        if (localz==(lattice::LZ - 1) && stencil::Ci_xyz(z)[Q]> 0) { //(note that the z direction goes from 0 to lattice::LZ-1)
                                                     //if the next lattice point in the z direction is divisible
                                                     //by lattice::LZ (so we are at z=lattice::LZ-1) and we are pointing in the +z
                                                     //direction

            neighbor += -(lattice::LZ - 1); //reduce k by lattice::LZ-1 so we are now at z=0

        }
        else if (localz==(0) && stencil::Ci_xyz(z)[Q] <0) { //if the current lattice point in the z direction is
                                                        //divisible by lattice::LZ (so we are at z=0) and we are pointing
                                                        //in the -z direction

            neighbor += (lattice::LZ - 1); //increase k by lattice::LZ-1 so we are now at z=lattice::LZ-1

        }
        else if (stencil::Ci_xyz(z)[Q] != 0) { //Else calculate neighbors normally

            neighbor += stencil::Ci_xyz(z)[Q]; //For z direction, the neighbor is just +Ci_z[Q]

        }
    }
    if(lattice::LY> 1) {
        int localY=computeY(lattice::LY,lattice::LZ,k);
        if (localY == (lattice::LY - 1) && stencil::Ci_xyz(y)[Q]> 0) { //(note that the y direction goes from 0 to lattice::LY-1)
                                                            //if the next lattice point in the y direction is
                                                            //divisible by lattice::LY (so we are at z=lattice::LY-1) and we are
                                                            //pointing in the +y direction

            neighbor += -(lattice::LZ) * (lattice::LY - 1);

        }
        else if (localY == 0 && stencil::Ci_xyz(y)[Q] <0) { //...

            neighbor += (lattice::LZ) * (lattice::LY - 1);

        }
        else if (stencil::Ci_xyz(y)[Q] != 0) {

            neighbor += stencil::Ci_xyz(y)[Q] * lattice::LZ;

        }
    }
    if(lattice::LXdiv> 1) {
        int localX=computeX(lattice::LY,lattice::LZ,k);
        if (localX == (lattice::LXdiv - 1) && stencil::Ci_xyz(x)[Q]> 0) { //...

            neighbor += -(lattice::LZ) * lattice::LY * (lattice::LXdiv - 1);

        }
        else if (localX == 0 && stencil::Ci_xyz(x)[Q] <0) {

            neighbor += (lattice::LZ) * lattice::LY * (lattice::LXdiv - 1);

        }
        else if (stencil::Ci_xyz(x)[Q] != 0) {

            neighbor += stencil::Ci_xyz(x)[Q] * lattice::LZ * lattice::LY;
            
        }
    }
    
    return k + neighbor; //return k + our neighbor offset
}

/**
 * \details This will iterate through the lattice and calculate the neighbors depending on whether the current
 *          lattice point is on a periodic boundary or not.
 */
template<class lattice, class stencil>
inline void Data_Base<lattice,stencil>::generateNeighbors() { //Loop over all lattice points and calculate the neghbor at each point

    #pragma omp parallel for schedule(guided)
    for (int k = 0; k <lattice::N; k++) { //For loop over all lattice points
        
        for(int q = 0; q < stencil::Q; q++) {

            if(!m_Geometry.isPeriodic(k)) { //If not periodic

                mv_Neighbors[k*stencil::Q + q] = getOneNeighbor(k, q);
                
            }
            else { //Else if periodic

                mv_Neighbors[k*stencil::Q + q] = getOneNeighborPeriodic(k, q);
                
            }

        }

    }
    
}

/**
 * \brief The Data1 class inherits from the Data_Base class but also provides information about the storage of
 *        distributions, streaming and communication of distributions.
 * Much of the functionality is the same as Data_Base, but we have an object of the distribution relevant to this
 * data type in the class. The stream() and getStreamIndex() function determine how streaming occurs for this data
 * type.
 * \tparam stencil Velocity Stencil of class using this data type.
 */
template<class lattice, class stencil>
class Data1 : public Data_Base<lattice, stencil> {
    private:

        /**
        * \brief This struct contains the distribution information relevant to this data class.
        * Distribution class will allocate memory to
        * distribution arrays and contains the
        * streamIndex function which is returns the
        * index of the neighboring lattice point in
        * the direction Q.
        */
        struct Distribution_Derived : public Distribution_Base<stencil> { //
            
            /**
             * \brief The constructor for the class.
             * This constructor will call the constructor for the base distribution class, calculate the opposite
             * indices at each index Q and allocate memory for the new and old distributions.
             * \param neighbors reference to a vector containing the neighboring lattice points at each point. Used to
             *                  construct the Distribution_Base class.
             */

            Distribution_Derived(std::vector<int>& neighbors) : Distribution_Base<stencil>(neighbors), mv_Neighbors(neighbors) { //Initialise mv_DistNeighbors

                Distribution_Base<stencil>::mv_Distribution.resize(stencil::Q * lattice::N); //Array size is number of
                                                                                  //directions times number of
                                                                                  //lattice points
                Distribution_Base<stencil>::mv_OldDistribution.resize(stencil::Q * lattice::N); //Old distributions needed
                                                                                     //in this case
                
            }

            Distribution_Derived(Distribution_Derived& other) : Distribution_Base<stencil>(other.mv_Neighbors), mv_Neighbors(other.mv_Neighbors) { //Initialise mv_DistNeighbors

                Distribution_Base<stencil>::mv_Distribution.resize(stencil::Q * lattice::N); //Array size is number of
                                                                                  //directions times number of
                                                                                  //lattice points
                Distribution_Base<stencil>::mv_OldDistribution.resize(stencil::Q * lattice::N); //Old distributions needed
                                                                                     //in this case
                
            }

            /**
             * \brief Returns the opposite index at the chosen index (Rotation by 180 degrees).
             */

            std::vector<int>& mv_Neighbors;

            /**
             * \brief Returns the index that the current distribution will be streamed to.
             * \details In this case, this just returns the neighbor at the current lattice point in the direction Q.
             * \param k Index of current lattice point.
             * \param Q Discrete velocity direction (e.g. 0-8 for D2Q9).
             * \return Index of distribution vector that the distribution will be streamed to.
             */
            inline int streamIndex(const int k, const int Q) {

                return Distribution_Base<stencil>::mv_DistNeighbors[k * stencil::Q + Q]; //Return neighbor of lattice point k in direction Q

            }

        };

        Distribution_Derived m_Distribution; //!<Object of distribution.
        
    public:
    
        /**
         * \brief Streaming step in the LBM algorithm.
         */
        inline void stream();

        /**
         * \brief This function streams the distributions to the neighboring processor.
         */
        inline void communicateDistribution();

        /**
         * \brief This constructor calls the constructor of the base disribution using the neighbor information.
         */
        Data1() : m_Distribution(Data_Base<lattice,stencil>::getInstance().mv_Neighbors) { //Construct distribution

        }

        Data1(Data1<lattice,stencil>& other) : m_Distribution(other.m_Distribution) { //Construct distribution

        }

        using DistributionData = Distribution_Derived; //!<Typedef so that the distribution class is available outside of this class.

        /**
         * \brief This function returns a reference to the distribution object stored within this class.
         * \return Reference to object of distribution.
         */
        inline Distribution_Derived& getDistributionObject() {

            return m_Distribution;

        }
        
};

/**
 * \details The stream() function does nothing in this class as streaming is implemented using the
 *          getStreamIndex() function in this data type.
 */
template<class lattice, class stencil>
inline void Data1<lattice, stencil>::stream() { //Not used in this data type

}

/**
 * \details This performs the communicateDistribution() function for the chosen parallelisation method, which
 *          should perform the streaming step across MPI boundaries.
 */
template<class lattice, class stencil>
inline void Data1<lattice, stencil>::communicateDistribution() {
    
    lattice::communicateDistribution(m_Distribution);
    
}
