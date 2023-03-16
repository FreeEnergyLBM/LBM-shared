#ifndef DATA_HEADER
#define DATA_HEADER
#include "Geometry.hh"
#include "Parameters.hh"
#include "Service.hh"

template<class stencil,int id>
class Data1{
    private:

        template<int i>
        struct Distribution_Derived:Distribution<stencil,i>{
            
            Distribution_Derived(std::vector<int>& neighbors):mv_DistNeighbors(neighbors){
                
                for(int idx=0;idx<stencil::Q;idx++){
                    opposites[idx]=stencil::Ci_xyz(x)[idx]*LZ*LY+stencil::Ci_xyz(y)[idx]*LZ+stencil::Ci_xyz(z)[idx];
                }

                Distribution<stencil,i>::mv_Distribution.resize(stencil::Q*N);
                Distribution<stencil,i>::mv_OldDistribution.resize(stencil::Q*N);
                
            }

            int opposites[stencil::Q];

            int streamIndex(const int k,const int Q) override{

                return mv_DistNeighbors[k*stencil::Q+Q];

            }

            std::vector<int>& mv_DistNeighbors;

        };

        int getOneNeighbor(const int k,const int Q);

        int getOneNeighborPeriodic(const int k,const int Q);

        int opp[stencil::Q];
        
        Geometry m_Geometry;

        std::vector<int> mv_Neighbors;

        enum{x=0,y=1,z=2};

        Distribution_Derived<id> m_Distribution;

    public:

        Data1():m_Distribution(mv_Neighbors){

            for(int idx=0;idx<stencil::Q;idx++){
                opp[idx]=stencil::Ci_xyz(x)[idx]*LZ*LY+stencil::Ci_xyz(y)[idx]*LZ+stencil::Ci_xyz(z)[idx];
            }

            mv_Neighbors.resize(stencil::Q*LX*LY*LZ);

            generateNeighbors();

        }

        void stream();

        void generateNeighbors();

        std::vector<int>& getNeighbors();

        Distribution_Derived<id>& getDistributionObject(){

            return m_Distribution;

        }

        int iterate(int k);

        int iterateFluid(int k);

};

template<typename stencil,int id>
void Data1<stencil,id>::stream(){

}

template<typename stencil,int id>
std::vector<int>& Data1<stencil,id>::getNeighbors(){
    
    return mv_Neighbors;

}

template<typename stencil,int id>
int Data1<stencil,id>::iterate(const int k){

    if (k>=N-1) return -1;
    else return k+1;

}

template<typename stencil,int id>
int Data1<stencil,id>::iterateFluid(const int k){

    if (k>=N-1) return -1;
    else if (m_Geometry.isSolid(k+1)) return iterateFluid(k+1);
    else return k+1;

}

template<typename stencil,int id>
int Data1<stencil,id>::getOneNeighbor(const int k,const int Q){
    
    return k+opp[Q];
        
}

template<typename stencil,int id>
int Data1<stencil,id>::getOneNeighborPeriodic(const int k,const int Q){

    int neighbor=0;

    if(LZ>1){
        if ((k+1)%(LZ)==0&&stencil::Ci_xyz(z)[Q]>0){

            neighbor+=-(LZ-1);

        }
        else if ((k)%(LZ)==0&&stencil::Ci_xyz(z)[Q]<0){

            neighbor+=(LZ-1);

        }
        else if (stencil::Ci_xyz(z)[Q]!=0){

            neighbor+=stencil::Ci_xyz(z)[Q];

        }
    }
    if(LY>1){
        if (((k)/(LZ)+1)%(LY)==0&&stencil::Ci_xyz(y)[Q]>0){

            neighbor+=-(LZ)*(LY-1);

        }
        else if (((k)/(LZ))%(LY)==0&&stencil::Ci_xyz(y)[Q]<0){

            neighbor+=(LZ)*(LY-1);

        }
        else if (stencil::Ci_xyz(y)[Q]!=0){

            neighbor+=stencil::Ci_xyz(y)[Q]*LZ;

        }
    }
    if(LX>1){
        if ((k/(LZ)/(LY)+1)%(LX)==0&&stencil::Ci_xyz(x)[Q]>0){

            neighbor+=-(LZ)*LY*(LX-1);

        }
        else if (((k)/(LZ)/(LY))%(LX)==0&&stencil::Ci_xyz(x)[Q]<0){

            neighbor+=(LZ)*LY*(LX-1);

        }
        else if (stencil::Ci_xyz(x)[Q]!=0){

            neighbor+=stencil::Ci_xyz(x)[Q]*LZ*LY;
            
        }
    }
    
    return k+neighbor;
}

template<typename stencil,int id>
void Data1<stencil,id>::generateNeighbors(){

    int k=0;

    while(k>=0){

        for(int q=0;q<stencil::Q;q++){

            if(!m_Geometry.isPeriodic(k)) {

                mv_Neighbors[k*stencil::Q+q]=getOneNeighbor(k,q);
                
            }
            else {

                mv_Neighbors[k*stencil::Q+q]=getOneNeighborPeriodic(k,q);

            }

        }

        k=iterate(k);
    }
}
#endif