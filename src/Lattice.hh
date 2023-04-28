#ifndef LATTICE_HEADER
#define LATTICE_HEADER

struct LatticeProperties{
    LatticeProperties(int lx,int ly, int lz=1, double DT=1.0):m_LX(lx),m_LY(ly),m_LZ(lz),m_N(lx*ly*lz),m_DT(DT),m_NDIM((lx<=1||ly<=1||lz<=1)*-1+3){
   
    }
    const int m_LX;
    int m_LXdiv;
    const int m_LY;
    const int m_LZ;
    int m_N;
    const int m_NDIM;
    const double m_DT;
    
};

#endif