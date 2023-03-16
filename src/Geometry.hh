#ifndef GEOMETRY_HEADER
#define GEOMETRY_HEADER
#include<map>
#include "Service.hh"

class Geometry{
    public:

        void initialise();

        bool isPeriodic(int k);

        bool isSolid(int k);

    private:

        enum{Periodic=3,Solid=2,Wall=1,Fluid=0};

        std::vector<int> geometry;


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
        LX==1) return true;
    
    else return false;
}

bool Geometry::isSolid(int k){

    int yAtCurrentk=computeY(k);

    if (yAtCurrentk<=1||yAtCurrentk>=LY-2) return true;
    
    else return false;

}

#endif