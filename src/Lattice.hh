#ifndef LATTICE_HEADER
#define LATTICE_HEADER

template<int lx, int ly,int lz=1>
struct LatticeProperties{
    LatticeProperties(double DT=1.0):m_N(lx*ly*lz),m_DT(DT){
   
    }
    static constexpr int m_LX=lx;
    int m_LXdiv;
    static constexpr int m_LY=ly;
    static constexpr int m_LZ=lz;
    int m_N;
    static constexpr int m_NDIM=(lx<=1||ly<=1||lz<=1)*-1+3;
    const double m_DT;
    
};

#endif