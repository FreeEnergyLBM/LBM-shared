#ifndef MPIPARALLEL_HEADER
#define MPIPARALLEL_HEADER
#ifdef MPIPARALLEL
#include <mpi.h>

template<int num_neighbors>
class Parallel{
    public:
        Parallel(){

            if(MAXNEIGHBORS<num_neighbors) MAXNEIGHBORS=num_neighbors;
        
            if (LX%NUMPROCESSORS==0) {
                LXdiv=(LX/NUMPROCESSORS+2*num_neighbors);
            }
            else{
                throw runtime_error(std::string("Currently, the number of cores must be divisible by the size of the domain in the x direction."));
            }
            /*
            else if (CURPROCESSOR<LX%NUMPROCESSORS){
                LXdiv=((LX-LX%NUMPROCESSORS)/NUMPROCESSORS+1+2*num_neighbors);
            }
            else {
                LXdiv=((LX-LX%NUMPROCESSORS)/NUMPROCESSORS+2*num_neighbors);
            }
            */
            N=LXdiv*LY*LZ;

        }
};

template<class stencil,int num_neighbors>
class X_Parallel:public Parallel<num_neighbors>{
    public:

        X_Parallel();

        template<class parameter>
        static void communicate(parameter& obj);

        template<class distribution>
        static void communicateDistribution(distribution& obj);

        

    private:

        static int m_LeftNeighbor;
        static int m_RightNeighbor;
        static char* m_MPIBuffer;
        static MPI_Datatype DistributionVector;
        

};

template<class stencil,int num_neighbors>
int X_Parallel<stencil,num_neighbors>::m_LeftNeighbor;

template<class stencil,int num_neighbors>
int X_Parallel<stencil,num_neighbors>::m_RightNeighbor;

template<class stencil,int num_neighbors>
char* X_Parallel<stencil,num_neighbors>::m_MPIBuffer;

template<class stencil,int num_neighbors>
MPI_Datatype X_Parallel<stencil,num_neighbors>::DistributionVector;

template<class stencil,int num_neighbors>
template<class parameter>
void X_Parallel<stencil,num_neighbors>::communicate(parameter& obj){
    
    MPI_Isend(&obj.getParameter()[N*parameter::m_Num-(num_neighbors+1)*LY*LZ],num_neighbors*LY*LZ*parameter::m_Num,mpi_get_type<typename parameter::ParamType>(),m_RightNeighbor,0,MPI_COMM_WORLD,&request);
    MPI_Isend(&obj.getParameter()[num_neighbors*LY*LZ*parameter::m_Num],num_neighbors*LY*LZ*parameter::m_Num,mpi_get_type<typename parameter::ParamType>(),m_LeftNeighbor,1,MPI_COMM_WORLD,&request);
    MPI_Irecv(&obj.getParameter()[N*parameter::m_Num-num_neighbors*LY*LZ*parameter::m_Num],num_neighbors*LY*LZ*parameter::m_Num,mpi_get_type<typename parameter::ParamType>(),m_RightNeighbor,1,MPI_COMM_WORLD,&request);
    MPI_Irecv(&obj.getParameter()[0],num_neighbors*LY*LZ*parameter::m_Num,mpi_get_type<typename parameter::ParamType>(),m_LeftNeighbor,0,MPI_COMM_WORLD,&request);
    MPI_Barrier(MPI_COMM_WORLD);
}

template<class stencil,int num_neighbors>
template<class distribution>
void X_Parallel<stencil,num_neighbors>::communicateDistribution(distribution& obj){

    if (CURPROCESSOR<NUMPROCESSORS){
        int leftorright;
        int id=0;
        for(int idx=1; idx<stencil::Q;idx++){
            leftorright=stencil::Ci_xyz(0)[idx];
            if(leftorright==-1){
            
                MPI_Isend(&obj.getDistribution()[(num_neighbors-1)*LY*LZ*stencil::Q+idx],1,DistributionVector,m_LeftNeighbor,id,MPI_COMM_WORLD,&request);
                
            }
            else if(leftorright==1){
            
                MPI_Isend(&obj.getDistribution()[N*stencil::Q-(num_neighbors)*LY*LZ*stencil::Q+idx],1,DistributionVector,m_LeftNeighbor,id,MPI_COMM_WORLD,&request);
                
            }
            id+=1;
        }
        id=0;
        for(int idx=1; idx<stencil::Q;idx++){
            leftorright=stencil::Ci_xyz(0)[idx];
            if(leftorright==-1){
            
                MPI_Irecv(&obj.getDistribution()[N*stencil::Q-(num_neighbors+1)*LY*LZ*stencil::Q+idx],1,DistributionVector,m_RightNeighbor,id,MPI_COMM_WORLD,&request);
                
            }
            else if(leftorright==1){
            
                MPI_Irecv(&obj.getDistribution()[(num_neighbors)*LY*LZ*stencil::Q+idx],1,DistributionVector,m_RightNeighbor,id,MPI_COMM_WORLD,&request);
                
            }
            id+=1;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    
}

template<class stencil,int num_neighbors>
X_Parallel<stencil,num_neighbors>::X_Parallel(){

    const int bufSize=(LY*LZ*num_neighbors*(5+2)*2*2+1000)*sizeof(double);
    
    if(bufSize>MPIBUFFERSIZE){
        if(MPIBUFFERSIZE!=0)MPI_Buffer_detach(MPIBUFFER,&MPIBUFFERSIZE);
        if(MPIBUFFERSIZE!=0)delete[] MPIBUFFER;
        MPIBUFFER=new char[bufSize];
        MPI_Buffer_attach(MPIBUFFER,bufSize);
        MPIBUFFERSIZE=bufSize;
    }
    //std::cout<<MPIBUFFERSIZE<<std::endl;
    m_LeftNeighbor=CURPROCESSOR-1;
    if (m_LeftNeighbor == -1) m_LeftNeighbor=NUMPROCESSORS-1;
    m_RightNeighbor=CURPROCESSOR+1;
    if (m_RightNeighbor == NUMPROCESSORS) m_RightNeighbor=0;

    MPI_Type_vector(LZ*LY,1,stencil::Q,mpi_get_type<double>(),&DistributionVector);
    MPI_Type_commit(&DistributionVector);
    
}
#else
class No_Parallel{

};
#endif
#endif