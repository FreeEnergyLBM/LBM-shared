#ifndef CXYZGRADIENT_HEADER
#define CXYZGRADIENT_HEADER

template<class stencil, template<typename givenstencil> class parallel>
struct CentralXYZ{

    Data_Base<stencil,parallel<stencil>> m_Data;

    template<class parameter>
    inline double computeFirstDerivative(const parameter& val,const int direciton,const int k);

    template<class parameter>
    inline double computeLaplacian(const parameter& val,const int k);

    Geometry m_Geometry;
    
};

template<class stencil, template<typename givenstencil> class parallel>
template<class parameter>
double CentralXYZ<stencil,parallel>::computeFirstDerivative(const parameter& val,const int direction,const int k){

    double temp=0;
    for (int idx=0;idx<stencil::Q;idx++){
        
        if((m_Geometry.isSolid(m_Data.getNeighbors()[k*stencil::Q+idx]))){
            temp+=1.0/stencil::Cs2*stencil::Weights[idx]*stencil::Ci_xyz(direction)[idx]*0.5*(val.getParameter(m_Data.getNeighbors()[k*stencil::Q+stencil::Opposites[idx]]));
        }
        else{
            temp+=1.0/stencil::Cs2*stencil::Weights[idx]*stencil::Ci_xyz(direction)[idx]*0.5*(val.getParameter(m_Data.getNeighbors()[k*stencil::Q+idx]));
        }
        
    }
    return temp;

}

template<class stencil, template<typename givenstencil> class parallel>
template<class parameter>
double CentralXYZ<stencil,parallel>::computeLaplacian(const parameter& val,const int k){

    double temp=0;
    for (int idx=1;idx<stencil::Q;idx++){
    
        if((!m_Geometry.isSolid(m_Data.getNeighbors()[k*stencil::Q+idx]))) temp+=1.0/stencil::Cs2*stencil::Weights[idx]*2*(val.getParameter(m_Data.getNeighbors()[k*stencil::Q+idx])-val.getParameter(k));
        else temp+=1.0/stencil::Cs2*stencil::Weights[idx]*2*(val.getParameter(m_Data.getNeighbors()[k*stencil::Q+stencil::Opposites[idx]])-val.getParameter(k));

    }
    return temp;
}

#endif