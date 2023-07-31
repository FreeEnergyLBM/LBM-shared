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
 * TLattice point in every direction but it might be faster to just recalculate every time this is needed.
 * "Neighbor" refers to the TLattice point adjacent to any given latice point in a chosen discrete direction
 * Periodic boundaries work by setting the neighbor of each TLattice point at the edges of the domain so that the
 * top of the domain connects to the bottom, the left connects to the right etc.
 */


/**
 * \brief The Data_Base class provides information for the layout of neighboring TLattice points and communication
          for non-distribution parameters.
 * This class takes a TStencil as a template argument, as the velocity discretisation information and weights is 
 * needed. The class has public functions for communication, generating neighbors based on the
 * TStencil and a function to get a reference to the vector of neighbor indices.
 * \tparam TStencil Velocity Stencil of class using this data type.
 */
template<class TLattice, class TStencil>
class Data_Base{
    public:
        static Data_Base<TLattice,TStencil>& getInstance() {
            static Data_Base<TLattice,TStencil> instance;
            return instance;
        }
    private:
        //ADD VIRTUAL TO THIS
        /**
         * \brief This function returns the neighbor at the current TLattice point k in the direction Q.
         * \param k Index of current TLattice point.
         * \param Q Discrete velocity direction (e.g. 0-8 for D2Q9).
         * \return Lattice index of neighboring TLattice point in chosen direction.
         */
        inline int getOneNeighbor(const int k, const int Q);

        /**
         * \brief This function returns the neighbor at the current TLattice point k in the direction Q given that
         *        the TLattice point lies on a periodic boundary.
         * \param k Index of current TLattice point.
         * \param Q Discrete velocity direction (e.g. 0-8 for D2Q9).
         * \return Lattice index of neighboring TLattice point in chosen direction given that the current point is on a periodic boundary.
         */
        inline int getOneNeighborPeriodic(const int k, const int Q);

        std::array<int,TStencil::Q> OppositeOffset; //!<Opposite TLattice point offset in each direction.
        
        Geometry<TLattice> mGeometry; //!<Class containing geometry information (for periodic boundaries).

        std::vector<int> mv_Neighbors; //!<Vector containing neighbor information.

        enum{ x = 0, y = 1, z = 2 }; //!<Indices corresponding to x, y, z.
        
        template<class, class>
        friend class DataOldNew; //DataOldNew can access private members of the base class (will need to add new data
                            //types here but I will probably change how this works).

        template<class, class>
        friend class DataOldNewEquilibrium;

        /**
         * \brief The constructor for the class.
         * This constructor will calculate the opposite points at each index Q,
         * allocate memory for neighbors and fill the array of neighbors.
         */
        Data_Base() //Construct distribution
        {

            for(int idx = 0; idx <TStencil::Q; idx++) {

                OppositeOffset[idx] = TStencil::Ci_xyz(x)[idx] * TLattice::LZ * TLattice::LY + TStencil::Ci_xyz(y)[idx] * TLattice::LZ + TStencil::Ci_xyz(z)[idx];

            }

            mv_Neighbors.resize(TStencil::Q * TLattice::N); //Allocate memory for neighbors array
            
            generateNeighbors(); //Fill neighbors array
            
        }
        
    public:

        Data_Base(Data_Base<TLattice, TStencil> const& other) = delete;
        void operator=(Data_Base<TLattice, TStencil> const& other) = delete;

        using Stencil = TStencil; //!<Typedef so type info can be accessed from outside the class.

        /**
         * \brief This function communicates a chosen TParameter (halo regions are exchanged with neighboring
         *        processors).
         * \param obj Object of chosen TParameter.
         * \tparam TParameter type of object to be communicated.
         */
        template<class TParameter>
        inline void communicate(TParameter& obj);

        /**
         * \brief Function to fill neighbor array with neighbor information.
         */
        inline void generateNeighbors();

        /**
         * \brief Returns the neighbor array
         */
        inline std::vector<int>& getNeighbors();
        inline int& getNeighbor(int k,int q){return mv_Neighbors[k*TStencil::Q+q];};
        const inline std::vector<int>& getNeighbors() const;

};

/**
 * \details The function calls the communicate(obj) function from the parallel template class chosen when the
 *          class is used. This will perform the necessary communications so gradients etc. can be calculated
 *          across parallel regions.
 */
template<class TLattice, class TStencil>
template<class TParameter>
inline void Data_Base<TLattice,TStencil>::communicate(TParameter& obj) { //Not used in this data type

    //static_assert(is_base_of_template<Parameter,TParameter>::value,"ERROR: The object passed to this function cannot be communicated.");
    TLattice::communicate(obj);

}

/**
 * \details The function returns a reference to a vector containing the neighboring TLattice point for every point
 *          in every direction in the TStencil.
 */
template<class TLattice, class TStencil>
inline std::vector<int>& Data_Base<TLattice, TStencil>::getNeighbors() {
    
    return mv_Neighbors;

}

template<class TLattice, class TStencil>
const inline std::vector<int>& Data_Base<TLattice, TStencil>::getNeighbors() const {
    
    return mv_Neighbors;

}

/**
 * \details The neighbor of the current TLattice point is calculated from the current TLattice point + the offset
 *          in each direction, which is precomputed in the constructor and stored in a vector.
 */
template<class TLattice, class TStencil>
inline int Data_Base<TLattice,TStencil>::getOneNeighbor(const int k, const int Q) {
    
    return k + OppositeOffset[Q]; //The neighbor is the TLattice point plus the opposite offset in direction Q
        
}

/**
 * \details If we are at the first TLattice point, some of the adjacent points will be on the complete opposite
 *          side of the TLattice so we must account for this. The function will work out if we are on the edge of
 *          of the simulation in the x, y and z directions and apply offsets in each case.
 */
template<class TLattice, class TStencil>
inline int Data_Base<TLattice,TStencil>::getOneNeighborPeriodic(const int k, const int Q) {

    int neighbor = 0;

    if(TLattice::LZ> 1) {
        int localz=computeZ(TLattice::LY,TLattice::LZ,k);
        if (localz==(TLattice::LZ - 1) && TStencil::Ci_xyz(z)[Q]> 0) { //(note that the z direction goes from 0 to TLattice::LZ-1)
                                                     //if the next TLattice point in the z direction is divisible
                                                     //by TLattice::LZ (so we are at z=TLattice::LZ-1) and we are pointing in the +z
                                                     //direction

            neighbor += -(TLattice::LZ - 1); //reduce k by TLattice::LZ-1 so we are now at z=0

        }
        else if (localz==(0) && TStencil::Ci_xyz(z)[Q] <0) { //if the current TLattice point in the z direction is
                                                        //divisible by TLattice::LZ (so we are at z=0) and we are pointing
                                                        //in the -z direction

            neighbor += (TLattice::LZ - 1); //increase k by TLattice::LZ-1 so we are now at z=TLattice::LZ-1

        }
        else if (TStencil::Ci_xyz(z)[Q] != 0) { //Else calculate neighbors normally

            neighbor += TStencil::Ci_xyz(z)[Q]; //For z direction, the neighbor is just +Ci_z[Q]

        }
    }
    if(TLattice::LY> 1) {
        int localY=computeY(TLattice::LY,TLattice::LZ,k);
        if (localY == (TLattice::LY - 1) && TStencil::Ci_xyz(y)[Q]> 0) { //(note that the y direction goes from 0 to TLattice::LY-1)
                                                            //if the next TLattice point in the y direction is
                                                            //divisible by TLattice::LY (so we are at z=TLattice::LY-1) and we are
                                                            //pointing in the +y direction

            neighbor += -(TLattice::LZ) * (TLattice::LY - 1);

        }
        else if (localY == 0 && TStencil::Ci_xyz(y)[Q] <0) { //...

            neighbor += (TLattice::LZ) * (TLattice::LY - 1);

        }
        else if (TStencil::Ci_xyz(y)[Q] != 0) {

            neighbor += TStencil::Ci_xyz(y)[Q] * TLattice::LZ;

        }
    }
    if(TLattice::LXdiv> 1) {
        int localX=computeX(TLattice::LY,TLattice::LZ,k);
        if (localX == (TLattice::LXdiv - 1) && TStencil::Ci_xyz(x)[Q]> 0) { //...

            neighbor += -(TLattice::LZ) * TLattice::LY * (TLattice::LXdiv - 1);

        }
        else if (localX == 0 && TStencil::Ci_xyz(x)[Q] <0) {

            neighbor += (TLattice::LZ) * TLattice::LY * (TLattice::LXdiv - 1);

        }
        else if (TStencil::Ci_xyz(x)[Q] != 0) {

            neighbor += TStencil::Ci_xyz(x)[Q] * TLattice::LZ * TLattice::LY;
            
        }
    }
    
    return k + neighbor; //return k + our neighbor offset
}

/**
 * \details This will iterate through the TLattice and calculate the neighbors depending on whether the current
 *          TLattice point is on a periodic boundary or not.
 */
template<class TLattice, class TStencil>
inline void Data_Base<TLattice,TStencil>::generateNeighbors() { //Loop over all TLattice points and calculate the neghbor at each point

    #pragma omp parallel for schedule(guided)
    for (int k = 0; k <TLattice::N; k++) { //For loop over all TLattice points
        
        for(int q = 0; q < TStencil::Q; q++) {

            if(!mGeometry.isPeriodic(k)) { //If not periodic

                mv_Neighbors[k*TStencil::Q + q] = getOneNeighbor(k, q);
                
            }
            else { //Else if periodic

                mv_Neighbors[k*TStencil::Q + q] = getOneNeighborPeriodic(k, q);
                
            }

        }

    }
    
}

/**
 * \brief The DataOldNew class inherits from the Data_Base class but also provides information about the storage of
 *        distributions, streaming and communication of distributions.
 * Much of the functionality is the same as Data_Base, but we have an object of the distribution relevant to this
 * data type in the class. The stream() and getStreamIndex() function determine how streaming occurs for this data
 * type.
 * \tparam TStencil Velocity Stencil of class using this data type.
 */
template<class TLattice, class TStencil>
class DataOldNew : public Data_Base<TLattice, TStencil> {
    private:

        /**
        * \brief This struct contains the distribution information relevant to this data class.
        * Distribution class will allocate memory to
        * distribution arrays and contains the
        * streamIndex function which is returns the
        * index of the neighboring TLattice point in
        * the direction Q.
        */
        struct Distribution_Derived : public Distribution_Base<TStencil> { //
            
            /**
             * \brief The constructor for the class.
             * This constructor will call the constructor for the base distribution class, calculate the opposite
             * indices at each index Q and allocate memory for the new and old distributions.
             * \param neighbors reference to a vector containing the neighboring TLattice points at each point. Used to
             *                  construct the Distribution_Base class.
             */

            static constexpr bool saveEquilibrium = false;

            Distribution_Derived(std::vector<int>& neighbors) : Distribution_Base<TStencil>(neighbors), mv_Neighbors(neighbors) { //Initialise mv_DistNeighbors

                Distribution_Base<TStencil>::mv_Distribution.resize(TStencil::Q * TLattice::N); //Array size is number of
                                                                                  //directions times number of
                                                                                  //TLattice points
                Distribution_Base<TStencil>::mv_OldDistribution.resize(TStencil::Q * TLattice::N); //Old distributions needed
                                                                                     //in this case
                
            }

            Distribution_Derived(Distribution_Derived& other) : Distribution_Base<TStencil>(other.mv_Neighbors), mv_Neighbors(other.mv_Neighbors) { //Initialise mv_DistNeighbors

                Distribution_Base<TStencil>::mv_Distribution.resize(TStencil::Q * TLattice::N); //Array size is number of
                                                                                  //directions times number of
                                                                                  //TLattice points
                Distribution_Base<TStencil>::mv_OldDistribution.resize(TStencil::Q * TLattice::N); //Old distributions needed
                                                                                     //in this case
                
            }

            /**
             * \brief Returns the opposite index at the chosen index (Rotation by 180 degrees).
             */

            std::vector<int>& mv_Neighbors;

            /**
             * \brief Returns the index that the current distribution will be streamed to.
             * \details In this case, this just returns the neighbor at the current TLattice point in the direction Q.
             * \param k Index of current TLattice point.
             * \param Q Discrete velocity direction (e.g. 0-8 for D2Q9).
             * \return Index of distribution vector that the distribution will be streamed to.
             */
            inline int streamIndex(const int k, const int Q) {

                return Distribution_Base<TStencil>::mv_DistNeighbors[k * TStencil::Q + Q]; //Return neighbor of TLattice point k in direction Q

            }

        };

        Distribution_Derived mDistribution; //!<Object of distribution.
        
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
        DataOldNew() : mDistribution(Data_Base<TLattice,TStencil>::getInstance().mv_Neighbors) { //Construct distribution

        }

        DataOldNew(DataOldNew<TLattice,TStencil>& other) : mDistribution(other.mDistribution) { //Construct distribution

        }

        using DistributionData = Distribution_Derived; //!<Typedef so that the distribution class is available outside of this class.

        /**
         * \brief This function returns a reference to the distribution object stored within this class.
         * \return Reference to object of distribution.
         */
        inline Distribution_Derived& getDistributionObject() {

            return mDistribution;

        }
        
};

/**
 * \details The stream() function does nothing in this class as streaming is implemented using the
 *          getStreamIndex() function in this data type.
 */
template<class TLattice, class TStencil>
inline void DataOldNew<TLattice, TStencil>::stream() { //Not used in this data type

}

/**
 * \details This performs the communicateDistribution() function for the chosen parallelisation method, which
 *          should perform the streaming step across MPI boundaries.
 */
template<class TLattice, class TStencil>
inline void DataOldNew<TLattice, TStencil>::communicateDistribution() {
    
    TLattice::communicateDistribution(mDistribution);
    
}

template<class TLattice, class TStencil>
class DataOldNewEquilibrium : public Data_Base<TLattice, TStencil> {
    private:

        /**
        * \brief This struct contains the distribution information relevant to this data class.
        * Distribution class will allocate memory to
        * distribution arrays and contains the
        * streamIndex function which is returns the
        * index of the neighboring TLattice point in
        * the direction Q.
        */
        struct Distribution_Derived : public Distribution_Base<TStencil> { //
            
            /**
             * \brief The constructor for the class.
             * This constructor will call the constructor for the base distribution class, calculate the opposite
             * indices at each index Q and allocate memory for the new and old distributions.
             * \param neighbors reference to a vector containing the neighboring TLattice points at each point. Used to
             *                  construct the Distribution_Base class.
             */

            static constexpr bool SaveEquilibrium = true;

            Distribution_Derived(std::vector<int>& neighbors) : Distribution_Base<TStencil>(neighbors), mv_Neighbors(neighbors) { //Initialise mv_DistNeighbors

                Distribution_Base<TStencil>::mv_Distribution.resize(TStencil::Q * TLattice::N); //Array size is number of
                                                                                  //directions times number of
                                                                                  //TLattice points
                Distribution_Base<TStencil>::mv_OldDistribution.resize(TStencil::Q * TLattice::N); //Old distributions needed
                                                                                     //in this case
                mv_EquilibriumDistribution.resize(TStencil::Q * TLattice::N); //Old distributions needed
                                                                                     //in this case
                
            }

            Distribution_Derived(Distribution_Derived& other) : Distribution_Base<TStencil>(other.mv_Neighbors), mv_Neighbors(other.mv_Neighbors) { //Initialise mv_DistNeighbors

                Distribution_Base<TStencil>::mv_Distribution.resize(TStencil::Q * TLattice::N); //Array size is number of
                                                                                  //directions times number of
                                                                                  //TLattice points
                Distribution_Base<TStencil>::mv_OldDistribution.resize(TStencil::Q * TLattice::N); //Old distributions needed
                                                                                     //in this case
                mv_EquilibriumDistribution.resize(TStencil::Q * TLattice::N); //Old distributions needed
                                                                                     //in this case
                
            }

            /**
             * \brief Returns the opposite index at the chosen index (Rotation by 180 degrees).
             */

            std::vector<int>& mv_Neighbors;

            std::vector<int> mv_EquilibriumDistribution;

            /**
             * \brief Returns the index that the current distribution will be streamed to.
             * \details In this case, this just returns the neighbor at the current TLattice point in the direction Q.
             * \param k Index of current TLattice point.
             * \param Q Discrete velocity direction (e.g. 0-8 for D2Q9).
             * \return Index of distribution vector that the distribution will be streamed to.
             */
            inline int streamIndex(const int k, const int Q) {

                return Distribution_Base<TStencil>::mv_DistNeighbors[k * TStencil::Q + Q]; //Return neighbor of TLattice point k in direction Q

            }

            inline void saveEquilibriums(const double* equilibrium, int k) {

                std::copy(equilibrium,equilibrium+TStencil::Q,&mv_EquilibriumDistribution[k * TStencil::Q]); //Return neighbor of TLattice point k in direction Q

            }

        };

        Distribution_Derived mDistribution; //!<Object of distribution.
        
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
        DataOldNewEquilibrium() : mDistribution(Data_Base<TLattice,TStencil>::getInstance().mv_Neighbors) { //Construct distribution

        }

        DataOldNewEquilibrium(DataOldNewEquilibrium<TLattice,TStencil>& other) : mDistribution(other.mDistribution) { //Construct distribution

        }

        using DistributionData = Distribution_Derived; //!<Typedef so that the distribution class is available outside of this class.

        /**
         * \brief This function returns a reference to the distribution object stored within this class.
         * \return Reference to object of distribution.
         */
        inline Distribution_Derived& getDistributionObject() {

            return mDistribution;

        }
        
};

/**
 * \details The stream() function does nothing in this class as streaming is implemented using the
 *          getStreamIndex() function in this data type.
 */
template<class TLattice, class TStencil>
inline void DataOldNewEquilibrium<TLattice, TStencil>::stream() { //Not used in this data type

}

/**
 * \details This performs the communicateDistribution() function for the chosen parallelisation method, which
 *          should perform the streaming step across MPI boundaries.
 */
template<class TLattice, class TStencil>
inline void DataOldNewEquilibrium<TLattice, TStencil>::communicateDistribution() {
    
    TLattice::communicateDistribution(mDistribution);
    
}
