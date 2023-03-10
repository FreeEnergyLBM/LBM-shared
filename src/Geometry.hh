#ifndef GEOMETRY_HEADER
#define GEOMETRY_HEADER
#include<map>
class Geometry{
    public:
        Geometry(){
            //initialise();
        }

        void initialise();

        bool isPeriodic(int k);

    private:
        enum{Periodic=3,Solid=2,Wall=1,Fluid=0};
        std::vector<int> geometry;

};

bool Geometry::isPeriodic(int k){
    if (k%(LZ-1)==0||(k-1)%(LZ-1)==0||k/(LZ-1)%(LY-1)==0||((k)/(LZ-1)-1)%(LY-1)==0||(k/(LZ-1)/(LY-1))%(LX-1)==0||((k-1)/(LZ-1)/(LY-1)-1)%(LX-1)==0) return true;
    else return false;
}
#endif