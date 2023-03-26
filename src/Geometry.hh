#ifndef GEOMETRY_HEADER
#define GEOMETRY_HEADER
#include<map>
#include "Service.hh"

//Geometry.hh: This is not fully implemented but will contain information relevant to boundary conditions

class Geometry{
    public:

        void initialise(); //Perform necessary initialisation

        bool isPeriodic(int k); //Returns true if the current lattice point lies on a periodic boundary

        bool isSolid(int k); //Returns true if the current lattice point lies on a solid boundary

    private:

        enum{Periodic=3,Solid=2,Wall=1,Fluid=0}; //IDs for each boundary type

        std::vector<int> geometry; //Vector containing geometry information (not used right now)


};

bool Geometry::isPeriodic(int k){

    if (((k%(LZ-1)==0||
         (k-1)%(LZ-1)==0)&&LZ>1)||
        ((k/(LZ)%(LY-1)==0||
         ((k)/(LZ)-1)%(LY-1)==0)&&LY>1)||
        (((k/(LZ)/(LY))%(LX-1)==0||
          ((k)/(LZ)/(LY)-1)%(LX-1)==0||
          ((k)/(LZ)/(LY)-1)-1<0||
          ((k)/(LZ)/(LY)-1)+1>LX-1)&&
         LX>1)||
        LZ==1||
        LY==1||
        LX==1) return true; //These conditions just check whether the lattice point is on the edge of the 
                            //simulation domain
    
    else return false;
}

bool Geometry::isSolid(int k){

    int yAtCurrentk=computeY(k);

    if (yAtCurrentk<=1||yAtCurrentk>=LY-2) return true; //Change this condition to control where the solid is
    
    else return false;

}

#endif