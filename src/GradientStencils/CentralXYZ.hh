#ifndef CXYZGRADIENT_HEADER
#define CXYZGRADIENT_HEADER

template<class data>
struct CentralXYZ{

    data m_Data;

    template<class parameter>
    inline double computeFirstDerivative(const parameter& val,const int direciton,const int k);

    template<class parameter>
    inline double computeLaplacian(const parameter& val,const int k);

    Geometry m_Geometry;
    
};

template<class data>
template<class parameter>
double CentralXYZ<data>::computeFirstDerivative(const parameter& val,const int direction,const int k){

    double temp=0;
    for (int idx=0;idx<data::Stencil::Q;idx++){
        
        if((m_Geometry.isSolid(m_Data.getNeighbors()[k*data::Stencil::Q+idx]))){
            temp+=1.0/data::Stencil::Cs2*data::Stencil::Weights[idx]*data::Stencil::Ci_xyz(direction)[idx]*0.5*(val.getParameter(m_Data.getNeighbors()[k*data::Stencil::Q+data::Stencil::Opposites[idx]]));
        }
        else{
            temp+=1.0/data::Stencil::Cs2*data::Stencil::Weights[idx]*data::Stencil::Ci_xyz(direction)[idx]*0.5*(val.getParameter(m_Data.getNeighbors()[k*data::Stencil::Q+idx]));
        }
        
    }
    return temp;

}

template<class data>
template<class parameter>
double CentralXYZ<data>::computeLaplacian(const parameter& val,const int k){

    double temp=0;
    for (int idx=0;idx<data::Stencil::Q;idx++){
    
        if((!m_Geometry.isSolid(m_Data.getNeighbors()[k*data::Stencil::Q+idx]))) temp+=1.0/data::Stencil::Cs2*data::Stencil::Weights[idx]*2*(val.getParameter(m_Data.getNeighbors()[k*data::Stencil::Q+idx])-val.getParameter(k));

    }
    return temp;
}

#endif