#ifndef CHEMFORCE_HEADER
#define CHEMFORCE_HEADER
#include "../Parameters.hh"
#include<iostream>

//ExternalForce.hh: Contains the force class for a constant applied body force in a given direction. This is
//unfinished (should be able to specify magnitude and direction).

class ChemicalForce{
    public:
        template<int lx, int ly,int lz>
        ChemicalForce(LatticeProperties<lx,ly,lz>& properties):m_ChemicalPotential(properties),m_GradOrderParameter(properties),m_LaplacianOrderParameter(properties),m_OrderParameter(properties),m_Density(properties),DT(properties.m_DT){}

        double compute(int xyz,int k) const; //Return force at lattice point k in direction xyz

        void precompute(int k); //Perform any neccessary computations before force is computed

        double computeDensitySource(int k) const; //Calculate any possible source/correction term for density

        double computeVelocitySource(int xyz,int k) const; //Calculate any possible source/correction term for
                                                           //velocity

        void postprocess(int k); //Perform any necessary postprocessing

    private:

        const double& DT;
        double m_A=0.00015;
        double m_Kappa=0.0003;

        ChemicalPotential m_ChemicalPotential;

        GradientOrderParameter m_GradOrderParameter;

        LaplacianOrderParameter m_LaplacianOrderParameter;

        OrderParameter m_OrderParameter;

        Density m_Density; //Density

};

double ChemicalForce::compute(int xyz,int k) const{

    return 0;//m_ChemicalPotential.getParameter(k)*m_GradOrderParameter.getParameterPointer(k)[xyz];

}

void ChemicalForce::precompute(int k){ //Not necessary

    m_ChemicalPotential.getParameter(k)=-m_A*m_OrderParameter.getParameter(k)+m_A*m_OrderParameter.getParameter(k)*m_OrderParameter.getParameter(k)*m_OrderParameter.getParameter(k)-m_Kappa*m_LaplacianOrderParameter.getParameter(k);

}

void ChemicalForce::postprocess(int k){ //Not necessary
    
}

double ChemicalForce::computeDensitySource(int k) const{ //Not necessary

    return 0.0;

}

double ChemicalForce::computeVelocitySource(int xyz,int k) const{ //Need to correct velocity

    return +compute(xyz,k)*DT/(2.0*m_Density.getParameter(k));
    
}

#endif