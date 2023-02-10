#ifndef FORCES_HEADER
#define FORCES_HEADER
#include "Parameters.hh"

class BodyForce{
    public:
        double compute(int xyz) const;
        void precompute();
        double computeDensitySource() const;
        double computeVelocitySource(int xyz) const;
    private:
        double magnitude=0.01;
        Density<double> m_Density;
};

double BodyForce::compute(int xyz) const{
    return (xyz==0)*magnitude*m_Density.getParameter(0);
}

void BodyForce::precompute(){
    
}

double BodyForce::computeDensitySource() const{
    return 0.0;
}

double BodyForce::computeVelocitySource(int xyz) const{
    return +compute(xyz)*DT/2.0;
}
#endif