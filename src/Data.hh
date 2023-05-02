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
 * needed. It also takes the MPI parallelisation method as a template argument as the communication may change
 * based on the data layout. The class has public functions for communication, generating neighbors based on the
 * stencil and a function to get a reference to the vector of neighbor indices.
 */
template<class stencil,class parallel>
class Data_Base{
    private:
        
        /**
         * \brief This function returns the neighbor at the current lattice point k in the direction Q.
         * \param k Index of current lattice point.
         * \param Q Discrete velocity direction (e.g. 0-8 for D2Q9).
         */
        int getOneNeighbor(const int k,const int Q);

        /**
         * \brief This function returns the neighbor at the current lattice point k in the direction Q given that
         *        the lattice point lies on a periodic boundary.
         * \param k Index of current lattice point.
         * \param Q Discrete velocity direction (e.g. 0-8 for D2Q9).
         */
        int getOneNeighborPeriodic(const int k,const int Q);

        int OppositeOffset[stencil::Q]; //!< Opposite lattice point offset in each direction.
        
        Geometry m_Geometry; //!< Class containing geometry information (for periodic boundaries).

        std::vector<int> mv_Neighbors; //!< Vector containing neighbor information.

        enum{x=0,y=1,z=2}; //!< Indices corresponding to x, y, z.
        
        parallel m_Parallel; //!< Object of parallelisation class.

        template<class stencil1,class parallel1>
        friend class Data1; //Data1 can access private members of the base class (will need to add new data
                            //types here but I will probably change how this works).
        
    public:

        using Stencil=stencil; //!< Typedef so type info can be accessed from outside the class.

        #ifdef MPIPARALLEL
        /**
         * \brief This function communicates a chosen parameter (halo regions are exchanged with neighboring
         *        processors).
         * \param obj Object of chosen parameter.
         */
        template<class parameter>
        void communicate(parameter obj);
        #endif

        /**
         * \brief The constructor for the class.
         * This constructor will call the constructor for the given parallel class, calculate the opposite points
         * at each index Q, allocate memory for neighbors and fill the array of neighbors.
         */
        template<class prop>
        Data_Base(prop& properties) //Construct distribution
         : m_Geometry(properties),
           m_Parallel(properties),
           LX(properties.m_LX),
           LY(properties.m_LY),
           LZ(properties.m_LZ),
           LXdiv(properties.m_LXdiv),
           N(properties.m_N)
        {
            for(int idx=0;idx<stencil::Q;idx++){
                OppositeOffset[idx]=stencil::Ci_xyz(x)[idx]*LZ*LY+stencil::Ci_xyz(y)[idx]*LZ+stencil::Ci_xyz(z)[idx];
            }

            mv_Neighbors.resize(stencil::Q*properties.m_N); //Allocate memory for neighbors array
            
            generateNeighbors(); //Fill neighbors array
            
        }

        const int& LX;
        const int& LY;
        const int& LZ;
        const int& LXdiv;
        const int& N;
        /**
         * \brief Function to fill neighbor array with neighbor information.
         */
        void generateNeighbors();

        /**
         * \brief Returns the neighbor array
         */
        std::vector<int>& getNeighbors();

};

#ifdef MPIPARALLEL
/**
 * \details The function calls the communicate(obj) function from the parallel template class chosen when the
 *          class is used. This will perform the necessary communications so gradients etc. can be calculated
 *          across parallel regions.
 */
template<class stencil,class parallel>
template<class parameter>
void Data_Base<stencil,parallel>::communicate(parameter obj){ //Not used in this data type

    m_Parallel.communicate(obj);

}
#endif

/**
 * \details The function returns a reference to a vector containing the neighboring lattice point for every point
 *          in every direction in the stencil.
 */
template<class stencil,class parallel>
std::vector<int>& Data_Base<stencil,parallel>::getNeighbors(){
    
    return mv_Neighbors;

}

/**
 * \details The neighbor of the current lattice point is calculated from the current lattice point + the offset
 *          in each direction, which is precomputed in the constructor and stored in a vector.
 */
template<class stencil,class parallel>
int Data_Base<stencil,parallel>::getOneNeighbor(const int k,const int Q){
    
    return k+OppositeOffset[Q]; //The neighbor is the lattice point plus the opposite offset in direction Q
        
}

/**
 * \details If we are at the first lattice point, some of the adjacent points will be on the complete opposite
 *          side of the lattice so we must account for this. The function will work out if we are on the edge of
 *          of the simulation in the x, y and z directions and apply offsets in each case.
 */
template<class stencil,class parallel>
int Data_Base<stencil,parallel>::getOneNeighborPeriodic(const int k,const int Q){ 

    int neighbor=0;

    if(LZ>1){
        if ((k+1)%(LZ)==0&&stencil::Ci_xyz(z)[Q]>0){ //(note that the z direction goes from 0 to LZ-1)
                                                     //if the next lattice point in the z direction is divisible
                                                     //by LZ (so we are at z=LZ-1) and we are pointing in the +z
                                                     //direction

            neighbor+=-(LZ-1); //reduce k by LZ-1 so we are now at z=0

        }
        else if ((k)%(LZ)==0&&stencil::Ci_xyz(z)[Q]<0){ //if the current lattice point in the z direction is
                                                        //divisible by LZ (so we are at z=0) and we are pointing
                                                        //in the -z direction

            neighbor+=(LZ-1); //increase k by LZ-1 so we are now at z=LZ-1

        }
        else if (stencil::Ci_xyz(z)[Q]!=0){ //Else calculate neighbors normally

            neighbor+=stencil::Ci_xyz(z)[Q]; //For z direction, the neighbor is just +Ci_z[Q]

        }
    }
    if(LY>1){
        if (((k)/(LZ)+1)%(LY)==0&&stencil::Ci_xyz(y)[Q]>0){ //(note that the y direction goes from 0 to LY-1)
                                                            //if the next lattice point in the y direction is
                                                            //divisible by LY (so we are at z=LY-1) and we are
                                                            //pointing in the +y direction

            neighbor+=-(LZ)*(LY-1);

        }
        else if (((k)/(LZ))%(LY)==0&&stencil::Ci_xyz(y)[Q]<0){ //...

            neighbor+=(LZ)*(LY-1);

        }
        else if (stencil::Ci_xyz(y)[Q]!=0){

            neighbor+=stencil::Ci_xyz(y)[Q]*LZ;

        }
    }
    if(LXdiv>1){
        if ((k/(LZ)/(LY)+1)%(LXdiv)==0&&stencil::Ci_xyz(x)[Q]>0){ //...

            neighbor+=-(LZ)*LY*(LXdiv-1);

        }
        else if (((k)/(LZ)/(LY))%(LXdiv)==0&&stencil::Ci_xyz(x)[Q]<0){

            neighbor+=(LZ)*LY*(LXdiv-1);

        }
        else if (stencil::Ci_xyz(x)[Q]!=0){

            neighbor+=stencil::Ci_xyz(x)[Q]*LZ*LY;
            
        }
    }
    
    return k+neighbor; //return k + our neighbor offset
}

/**
 * \details This will iterate through the lattice and calculate the neighbors depending on whether the current
 *          lattice point is on a periodic boundary or not.
 */
template<class stencil,class parallel>
void Data_Base<stencil,parallel>::generateNeighbors(){ //Loop over all lattice points and calculate the neghbor at each point

    #ifdef OMPPARALLEL
    #pragma omp parallel for schedule( dynamic )
    #endif
    for (int k=0;k<N;k++){ //For loop over all lattice points
        
        for(int q=0;q<stencil::Q;q++){

            if(!m_Geometry.isPeriodic(k)) { //If not periodic

                mv_Neighbors[k*stencil::Q+q]=getOneNeighbor(k,q);
                
            }
            else { //Else if periodic

                mv_Neighbors[k*stencil::Q+q]=getOneNeighborPeriodic(k,q);

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
 */
template<class stencil,class parallel>
class Data1:public Data_Base<stencil,parallel>{
    private:

        /**
        * \brief This struct contains the distribution information relevant to this data class.
        * Distribution class will allocate memory to
        * distribution arrays and contains the
        * streamIndex function which is returns the
        * index of the neighboring lattice point in
        * the direction Q.
        */
        struct Distribution_Derived:public Distribution_Base<stencil>{ //
            
            /**
             * \brief The constructor for the class.
             * This constructor will call the constructor for the base distribution class, calculate the opposite
             * indices at each index Q and allocate memory for the new and old distributions.
             * \param neighbors reference to a vector containing the neighboring lattice points at each point. Used to
             *                  construct the Distribution_Base class.
             */
            template<class prop>
            Distribution_Derived(prop& properties,std::vector<int>& neighbors):Distribution_Base<stencil>(neighbors){ //Initialise mv_DistNeighbors
                
                for(int idx=0;idx<stencil::Q;idx++){ //Calculate the k offset for the neighbors in each direction
                    ma_Opposites[idx]=stencil::Opposites[idx];
                }

                Distribution_Base<stencil>::mv_Distribution.resize(stencil::Q*properties.m_N); //Array size is number of
                                                                                  //directions times number of
                                                                                  //lattice points
                Distribution_Base<stencil>::mv_OldDistribution.resize(stencil::Q*properties.m_N); //Old distributions needed
                                                                                     //in this case
                
            }

            /**
             * \brief Returns the opposite index at the chosen index (Rotation by 180 degrees).
             */
            int getOpposite(int idx){

                return ma_Opposites[idx];

            }


            int ma_Opposites[stencil::Q]; //!< Array containing the opposite indices at each index (rotated by 180 degrees).

            /**
             * \brief Returns the index that the current distribution will be streamed to.
             * \details In this case, this just returns the neighbor at the current lattice point in the direction Q.
             * \param k Index of current lattice point.
             * \param Q Discrete velocity direction (e.g. 0-8 for D2Q9).
             */
            int streamIndex(const int k,const int Q){

                return Distribution_Base<stencil>::mv_DistNeighbors[k*stencil::Q+Q]; //Return neighbor of lattice point k in direction Q

            }

        };

        Distribution_Derived m_Distribution; //!< Object of distribution.
        
    public:
    
        /**
         * \brief Streaming step in the LBM algorithm.
         */
        void stream();

        #ifdef MPIPARALLEL
        /**
         * \brief This function streams the distributions to the neighboring processor.
         */
        void communicateDistribution();
        #endif

        /**
         * \brief This constructor calls the constructor of the base disribution using the neighbor information.
         */
        template<class prop>
        Data1(prop& properties):Data_Base<stencil,parallel>(properties),m_Distribution(properties,Data_Base<stencil,parallel>::mv_Neighbors){ //Construct distribution

        }

        using DistributionData=Distribution_Derived; //!< Typedef so that the distribution class is available outside of this class.

        /**
         * \brief This function returns a reference to the distribution object stored within this class.
         */
        Distribution_Derived& getDistributionObject(){

            return m_Distribution;

        }
        
};

/**
 * \details The stream() function does nothing in this class as streaming is implemented using the
 *          getStreamIndex() function in this data type.
 */
template<class stencil,class parallel>
void Data1<stencil,parallel>::stream(){ //Not used in this data type

}

#ifdef MPIPARALLEL
/**
 * \details This performs the communicateDistribution() function for the chosen parallelisation method, which
 *          should perform the streaming step across MPI boundaries.
 */
template<class stencil,class parallel>
void Data1<stencil,parallel>::communicateDistribution(){
    #ifdef MPIPARALLEL
    Data_Base<stencil,parallel>::m_Parallel.communicateDistribution(m_Distribution);
    #endif
}
#endif
