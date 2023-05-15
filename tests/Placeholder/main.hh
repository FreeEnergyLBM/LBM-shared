#include<../../src/lbm.hh>

template<class lattice>
void setFluid(){
    OrderParameter<lattice> m_OrderParameter;
    for(int k = lattice::m_HaloSize; k < lattice::m_N-lattice::m_HaloSize; k++) {
        int yy=computeY(lattice::m_LY,lattice::m_LZ,k);
        if ((yy)>lattice::m_LY/2) m_OrderParameter.initialiseUser(-1.0,k);
    }
}

template<class lattice>
void setGeometry(){
    SolidLabels<lattice> m_SolidLabels;
    for(int k = lattice::m_HaloSize; k < lattice::m_N-lattice::m_HaloSize; k++) {
        int yAtCurrentk = computeY(lattice::m_LY, lattice::m_LZ, k);
        if (yAtCurrentk <= 1 || yAtCurrentk >= lattice::m_LY - 2 ) m_SolidLabels.initialiseUser(true,k);
    }
}