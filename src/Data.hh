#ifndef DATA_HEADER
#define DATA_HEADER
#include "Geometry.hh"
#include "Parameters.hh"
#include "Service.hh"

//Data.hh: Contains data class that will control how how data is accessed and how stremaing happens. These
//classes will contain derived distribution classes that will determine memory allocation for the distribution
//arrays and the streaming indices. The Data class also contains a vector of neighbors for each lattice point in
//every direction but it might be faster to just recalculate every time this is needed.
//
//"Neighbor" refers to the lattice point adjacent to any given latice point in a chosen discrete direction
//
//Periodic boundaries work by setting the neighbor of each lattice point at the edges of the domain so that the
//top of the domain connects to the bottom, the left connects to the right etc.

template<class stencil,class parallel> //Stencil information is needed as streaming indices and neighbors are determined by the
                        //velocity vectors
class Data_Base{
    private:
        

        int getOneNeighbor(const int k,const int Q); //Return Neighbor of lattice point k in direction Q

        int getOneNeighborPeriodic(const int k,const int Q); //Return Neighbor of lattice point k in direction Q
                                                             //given that k lieas on a periodic boundary

        int OppositeOffset[stencil::Q]; //Opposite lattice point offset in each direction
        
        Geometry m_Geometry; //Class containing geometry information (for periodic boundaries)

        std::vector<int> mv_Neighbors; //Vector containing neighbor information

        enum{x=0,y=1,z=2}; //Indices corresponding to x, y, z
        parallel m_Parallel;

        template<class stencil1,class parallel1>
        friend class Data1;
        
    public:
        using Stencil=stencil;
        template<class parameter>
        void communicate(parameter obj);

        template<class parameter>
        void initialiseHalos(parameter obj);

        Data_Base():m_Parallel(){ //Construct distribution

            for(int idx=0;idx<stencil::Q;idx++){
                OppositeOffset[idx]=stencil::Ci_xyz(x)[idx]*LZ*LY+stencil::Ci_xyz(y)[idx]*LZ+stencil::Ci_xyz(z)[idx];
            }

            mv_Neighbors.resize(stencil::Q*N); //Allocate memory for neighbors array

            generateNeighbors(); //Fill neighbors array

        }

        void stream(); //Optional seperate streaming function

        void generateNeighbors(); //Function to fill neighbor array with neighbor information

        std::vector<int>& getNeighbors(); //Returns the neighbor array

        int iterate(int k,bool all); //Increment k THIS WILL BE CHANGED

        int iterateFluid(int k, bool all); //Increment k only over fluid THIS WILL BE CHANGED

        int iterateFluid0(int k, bool all); //Increment k only over fluid THIS WILL BE CHANGED


        int iterateSolid(int k, bool all); //Increment k only over solid THIS WILL BE CHANGED

        int iterateSolid0(int k, bool all); //Increment k only over solid THIS WILL BE CHANGED


};

template<class stencil,class parallel>
template<class parameter>
void Data_Base<stencil,parallel>::communicate(parameter obj){ //Not used in this data type

    parallel::communicate(obj);

}

template<class stencil,class parallel>
void Data_Base<stencil,parallel>::stream(){ //Not used in this data type

}

template<class stencil,class parallel>
std::vector<int>& Data_Base<stencil,parallel>::getNeighbors(){
    
    return mv_Neighbors; //Return the vector containing neighbor information.

}

template<class stencil,class parallel>
int Data_Base<stencil,parallel>::iterate(const int k,bool all){

    if (k>=N-MAXNEIGHBORS*LY*LZ*(!all)-1) return -1; //If k is at the final lattice point, return -1 which will terminate the loop
    else return k+1; //Else return k incremented by one

}

template<class stencil,class parallel>
int Data_Base<stencil,parallel>::iterateFluid(const int k,bool all){

    if (k>=N-MAXNEIGHBORS*LY*LZ*(!all)-1) return -1; //If k is at the final lattice point, return -1 which will terminate the loop
    else if (m_Geometry.isSolid(k+1)) return iterateFluid(k+1,all); //Else if we are on a solid node, skip this
                                                                //and try again for the next node
    else return k+1; //Else return k incremented by one

}

template<class stencil,class parallel>
int Data_Base<stencil,parallel>::iterateFluid0(const int k,bool all){

    if (k>=N-MAXNEIGHBORS*LY*LZ*(!all)-1) return -1; //If k is at the final lattice point, return -1 which will terminate the loop
    else if (m_Geometry.isSolid(k+1)) return iterateFluid(k+1,all); //Else if we are on a solid node, skip this
                                                                //and try again for the next node
    else return k; //Else return k incremented by one

}

template<class stencil,class parallel>
int Data_Base<stencil,parallel>::iterateSolid(const int k,bool all){

    if (k>=N-MAXNEIGHBORS*LY*LZ*(!all)-1) return -1; //If k is at the final lattice point, return -1 which will terminate the loop
    else if (!m_Geometry.isSolid(k+1)) return iterateSolid(k+1,all); //Else if we are on a solid node, skip this
                                                                //and try again for the next node
    else return k+1; //Else return k incremented by one

}

template<class stencil,class parallel>
int Data_Base<stencil,parallel>::iterateSolid0(const int k,bool all){

    if (k>=N-MAXNEIGHBORS*LY*LZ*(!all)-1) return -1; //If k is at the final lattice point, return -1 which will terminate the loop
    else if (!m_Geometry.isSolid(k+1)) return iterateSolid(k+1,all); //Else if we are on a solid node, skip this
                                                                //and try again for the next node
    else return k; //Else return k incremented by one

}

template<class stencil,class parallel>
int Data_Base<stencil,parallel>::getOneNeighbor(const int k,const int Q){
    
    return k+OppositeOffset[Q]; //The neighbor is the lattice point plus the opposite offset in direction Q
        
}

template<class stencil,class parallel>
int Data_Base<stencil,parallel>::getOneNeighborPeriodic(const int k,const int Q){ //This function will calculate the neighbors
                                                                     //given that "k" lies on a periodic boundary.
                                                                     //For instance, if we are at the first
                                                                     //lattice point, some of the adjacent points
                                                                     //will be on the complete opposite side of
                                                                     //the lattice so we must account for this.

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

template<class stencil,class parallel>
void Data_Base<stencil,parallel>::generateNeighbors(){ //Loop over all lattice points and calculate the neghbor at each point

    int k=0;

    while(k>=0){ //While look over all lattice points

        for(int q=0;q<stencil::Q;q++){

            if(!m_Geometry.isPeriodic(k)) { //If not periodic

                mv_Neighbors[k*stencil::Q+q]=getOneNeighbor(k,q);
                
            }
            else { //Else if periodic

                mv_Neighbors[k*stencil::Q+q]=getOneNeighborPeriodic(k,q);

            }

        }

        k=iterate(k,true); //Increment k
    }
}

template<class stencil,class parallel> //Stencil information is needed as streaming indices and neighbors are determined by the
                        //velocity vectors
class Data1:public Data_Base<stencil,parallel>{
    private:

        struct Distribution_Derived:public Distribution_Base<stencil>{ //Distribution class will allocate memory to
                                                                //distribution arrays and contains the
                                                                //streamIndex function which is returns the
                                                                //index of the neighboring lattice point in
                                                                //the direction Q.
            
            Distribution_Derived(std::vector<int>& neighbors):mv_DistNeighbors(neighbors){ //Initialise mv_DistNeighbors
                
                for(int idx=0;idx<stencil::Q;idx++){ //Calculate the k offset for the neighbors in each direction
                    ma_Opposites[idx]=stencil::Opposites[idx];
                }

                Distribution_Base<stencil>::mv_Distribution.resize(stencil::Q*N); //Array size is number of
                                                                                  //directions times number of
                                                                                  //lattice points
                Distribution_Base<stencil>::mv_OldDistribution.resize(stencil::Q*N); //Old distributions needed
                                                                                     //in this case
                
            }

            int getOpposite(int idx) override{

                return ma_Opposites[idx];

            }


            int ma_Opposites[stencil::Q];

            int streamIndex(const int k,const int Q) override{

                return mv_DistNeighbors[k*stencil::Q+Q]; //Return neighbor of lattice point k in direction Q

            }

            std::vector<int>& mv_DistNeighbors; //Reference to vector containing neighbor information

        };

        Distribution_Derived m_Distribution; //Object of distribution
        
    public:

        void communicateDistribution();

        Data1():m_Distribution(Data_Base<stencil,parallel>::mv_Neighbors){ //Construct distribution

        }

        Distribution_Derived& getDistributionObject(){

            return m_Distribution; //Returns the distribution object stored in the class

        }
        
};

template<class stencil,class parallel>
void Data1<stencil,parallel>::communicateDistribution(){ //Not used in this data type
    
    parallel::communicateDistribution(m_Distribution);
    
}
#endif