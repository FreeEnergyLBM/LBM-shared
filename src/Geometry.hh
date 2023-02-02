
class Geometry{
    public:
        Geometry(){
            geometry[0]=Fluid;
        }
    private:
        enum{Solid=2,Wall=1,Fluid=0};
        int geometry[1];
};