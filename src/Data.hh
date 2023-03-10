#ifndef DATA_HEADER
#define DATA_HEADER
#include "Geometry.hh"
#include "Parameters.hh"
template<class stencil>
class Data1{
    private:
        template<int i>
        struct Distribution_Derived:Distribution<stencil,i>{
            Distribution_Derived(){
                for(int idx=0;idx<stencil::Q;idx++){
                    opposites[idx]=stencil::Ci_xyz(x)[idx]*LZ*LY+stencil::Ci_xyz(y)[idx]*LZ+stencil::Ci_xyz(z)[idx];
                }
                Distribution<stencil,i>::mv_Distribution.resize(stencil::Q*N);
                Distribution<stencil,i>::mv_OldDistribution.resize(stencil::Q*N);
            }
            int opposites[stencil::Q];
            int streamIndex(const int k,const int Q) override;

        };

        int getOneNeighbor(const int k,const int Q);
        int getOneNeighborPeriodic(const int k,const int Q);
        int opp[stencil::Q];
        Geometry m_Geometry;
        std::vector<int> mv_Neighbors;
        enum{x=0,y=1,z=2};
        Distribution_Derived<1> m_Distribution;
    public:
        Data1(){
            for(int idx=0;idx<stencil::Q;idx++){
                opp[idx]=stencil::Ci_xyz(x)[idx]*LZ*LY+stencil::Ci_xyz(y)[idx]*LZ+stencil::Ci_xyz(z)[idx];
            }
            mv_Neighbors.resize(stencil::Q*LX*LY*LZ);
            generateNeighbors();
        }
        void stream();

        void generateNeighbors();

        void getNeighbor(const int k);

        template<int i>
        Distribution_Derived<i> getDistributionObject(){
            return m_Distribution;
        }

        int iterate(int k);

};

template<typename stencil>
void Data1<stencil>::stream(){

}

template<typename stencil>
int Data1<stencil>::iterate(const int k){
    if (k>=N-1) return -1;
    else return k+1;
}

template<typename stencil>
int Data1<stencil>::getOneNeighbor(const int k,const int Q){
    
    return k+opp[Q];
        
}

template<typename stencil>
int Data1<stencil>::getOneNeighborPeriodic(const int k,const int Q){
        int neighbor=0;
        if(LZ>1){
            if (k%(LZ-1)==0&&stencil::Ci_xyz(z)[Q]>0){
                neighbor+=-(LZ-1);
            }
            else if ((k-1)%(LZ-1)==0&&stencil::Ci_xyz(z)[Q]<0){
                neighbor+=(LZ-1);
            }
        }
        if(LY>1){
            if (k/(LZ-1)%(LY-1)==0&&stencil::Ci_xyz(y)[Q]>0){
                neighbor+=-(LZ)*(LY-1);
            }
            else if (((k)/(LZ-1)-1)%(LY-1)==0&&stencil::Ci_xyz(y)[Q]<0){
                neighbor+=(LZ)*(LY-1);
            }
        }
        if(LX>1){
            if ((k/(LZ-1)/(LY-1))%(LX-1)==0&&stencil::Ci_xyz(x)[Q]>0){
                neighbor+=-(LZ)*LY*(LX-1);
            }
            else if (((k-1)/(LZ-1)/(LY-1)-1)%(LX-1)==0&&stencil::Ci_xyz(x)[Q]<0){
                neighbor+=(LZ)*LY*(LX-1);
            }
        }
        return neighbor+getOneNeighbor(k,Q);
}

template<typename stencil>
void Data1<stencil>::generateNeighbors(){
    int k=0;
    while(k>0){
        for(int q=0;q<stencil::Q;q++){
            if(!m_Geometry.isPeriodic(k)) mv_Neighbors[k]=getOneNeighbor(k,q);
            else mv_Neighbors[k]=getOneNeighbor(k,q);

            k=iterate(k);
        }
    }
}

template<typename stencil>
template<int i>
int Data1<stencil>::template Distribution_Derived<i>::streamIndex(const int k,const int Q){
    return k+opposites[Q];
}
#endif