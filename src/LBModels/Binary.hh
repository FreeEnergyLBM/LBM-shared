#pragma once
#include "../Collide.hh"
#include "../Parameters.hh"
#include "../Data.hh"
#include "../Parallel.hh"
#include "ModelBase.hh"
#include <utility>

//Binary.hh: Contains the details of the LBM model to solve an equation for phase separation. Each
//Model is given a "traits" class that contains stencil, data, force and boundary information

template<class properties>
struct DefaultTraitBinary{

    using Stencil = std::conditional_t<std::remove_reference<properties>::type::m_NDIM == 2, D2Q9, D3Q19>; //Here, D refers to the number of cartesian dimensions

    using Boundaries = std::tuple<BounceBack>;

    using Forces = std::tuple<OrderParameterGradients<CentralXYZ<properties, Stencil>>>;

    using Properties = properties;

};


template<class traits = DefaultTraitBinary<decltype(GETPROPERTIES())>>
class Binary: public CollisionBase<typename traits::Stencil>, public ModelBase<traits> { //Inherit from base class to avoid repetition of common
                                                      //calculations
                                                      
    public:

        void collide() override; //Collision step

        void initialise() override; //Initialisation step

        void computeMomenta() override; //Momenta (density, velocity) calculation

        const double& getDensity(const int k) const; //Return density at lattice point k

        const std::vector<double>& getVelocity() const; //Return vector of velocity

    private:

        double computeEquilibrium(const double& orderparam, const double* velocity,
                                   const int idx, const int k) const; //Calculate equilibrium in direction idx with a given
                                                        //density and velocity

        double computeModelForce(int xyz, int k) const; //Calculate forces specific to the model in direction xyz

        double computeCollisionQ(double& sum, int k, const double& old, const double& orderparam,
                                  const double* velocity, const int idx) const; //Calculate collision
                                                                                           //at index idx

        double computeOrderParameter(const double* distribution, const int k) const; //Calculate the order parameter
                                                                              //corresponding to the relative
                                                                              //concentrations of each phase

        static constexpr double m_Tau = 1.0; //TEMPORARY relaxation time
        static constexpr double m_InverseTau = 1.0 / m_Tau; //TEMPORARY inverse relaxation time

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions
        
        double m_Gamma = 1;

        OrderParameter m_OrderParameter; //Order Parameter
        Velocity m_Velocity; //Velocity
        ChemicalPotential m_ChemicalPotential;
        GradientOrderParameter m_GradOrderParameter;
        LaplacianOrderParameter m_LaplacianOrderParameter;

        vector<double>& orderparameter = m_OrderParameter.getParameter(); //Reference to vector of order parameters
        vector<double>& velocity = m_Velocity.getParameter(); //Reference to vector of velocities

};

template<class traits>
const double& Binary<traits>::getDensity(const int k) const { //This needs to be renamed

    return orderparameter[k]; //Return reference to order parameter at point k

}

template<class traits>
const std::vector<double>& Binary<traits>::getVelocity() const {

    return velocity; //Return reference to velocity vector

}

template<class traits>
void Binary<traits>::collide() {

    #pragma omp parallel for schedule(guided)
    for (int k = GETPROPERTIES().m_HaloSize; k <GETPROPERTIES().m_N - GETPROPERTIES().m_HaloSize; k++){ //loop over k

        double* old_distribution = ModelBase<traits>::m_Distribution.getDistributionOldPointer(k);

        double equilibriumsum = 0;

        for (int idx = traits::Stencil::Q - 1; idx>= 0; --idx) { //loop over discrete velocity directions
            //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
            //"computeCollisionQ"

            double collision = computeCollisionQ(equilibriumsum, k, old_distribution[idx], orderparameter[k], &velocity[k * traits::Stencil::D], idx);

            ModelBase<traits>::m_Distribution.getDistributionPointer(ModelBase<traits>::m_Distribution.streamIndex(k, idx))[idx] = collision;

        }
        
    }

    #ifdef MPIPARALLEL
    ModelBase<traits>::m_Data.communicateDistribution();
    #endif

}

template<class traits>
void Binary<traits>::initialise() { //Initialise model

    ModelBase<traits>::m_Data.generateNeighbors(); //Fill array of neighbor values (See Data.hh)
    
    #pragma omp parallel for schedule(guided) 
    for (int k=GETPROPERTIES().m_HaloSize; k<GETPROPERTIES().m_N - GETPROPERTIES().m_HaloSize; k++) { //loop over k

        double* distribution = ModelBase<traits>::m_Distribution.getDistributionPointer(k);
        double* old_distribution=ModelBase<traits>::m_Distribution.getDistributionOldPointer(k);

        m_ChemicalPotential.getParameter(k) = 0;

        int xx = computeX(GETPROPERTIES().m_LY, GETPROPERTIES().m_LZ, k);

        if (xx>=GETPROPERTIES().m_LX/2) orderparameter[k] = 1.0; //Set order parameter to 1 initially (This will change)
        else orderparameter[k] = -1.0;

        double equilibriumsum = 0;

        for (int idx = traits::Stencil::Q - 1; idx>= 0; idx--) {

            double equilibrium;

            if (idx> 0) equilibrium = computeEquilibrium(orderparameter[k], &velocity[k * traits::Stencil::D], idx, k);
            else equilibrium = orderparameter[k] - equilibriumsum;

            distribution[idx] = equilibrium; //Set distributions to equillibrium
            old_distribution[idx] = equilibrium;        

        }

    }

}


template<class traits>
void Binary<traits>::computeMomenta() { //Calculate order parameter

    #pragma omp parallel for schedule(guided)
    for (int k = GETPROPERTIES().m_HaloSize; k <GETPROPERTIES().m_N - GETPROPERTIES().m_HaloSize; k++) { //Loop over k

        double* distribution = ModelBase<traits>::m_Distribution.getDistributionPointer(k);

        orderparameter[k] = computeOrderParameter(distribution, k);

    }

    #ifdef MPIPARALLEL
    ModelBase<traits>::m_Data.communicate(m_OrderParameter);
    #endif
    
}

template<class traits>
double Binary<traits>::computeCollisionQ(double& equilibriumsum, const int k, const double& old, const double& orderparam,
                                            const double* velocity, const int idx) const {
                                        //Calculate collision step at a given velocity index at point k
    
    double forcexyz[traits::Stencil::D]; //Temporary array storing force in each cartesian direction

    //Force is the sum of model forces and given forces
    for(int xyz = 0; xyz <traits::Stencil::D; xyz++) forcexyz[xyz] = computeModelForce(xyz, k) + ModelBase<traits>::computeForces(xyz, k);
    
    //Sum of collision + force contributions
    if (idx> 0) {

        double eq = CollisionBase<typename traits::Stencil>::collideSRT(old, computeEquilibrium(orderparam, velocity, idx, k), m_InverseTau)
                    + CollisionBase<typename traits::Stencil>::forceGuoSRT(forcexyz, velocity, m_InverseTau, idx);

        equilibriumsum += eq;
        
        return eq;
    }
    else return m_OrderParameter.getParameter(k) - equilibriumsum + CollisionBase<typename traits::Stencil>::forceGuoSRT(forcexyz, velocity, m_InverseTau, idx);

}


template<class traits>
double Binary<traits>::computeEquilibrium(const double& orderparam, const double* velocity, const int idx, const int k) const {

    return traits::Stencil::Weights[idx] * (m_ChemicalPotential.getParameter(k) * m_Gamma / traits::Stencil::Cs2 + orderparam * CollisionBase<typename traits::Stencil>::computeVelocityFactor(velocity, idx));

}

template<class traits>
double Binary<traits>::computeModelForce(int k,int xyz) const {
    
    return 0;

}


template<class traits>
double Binary<traits>::computeOrderParameter(const double* distribution, const int k) const {//Order parameter calculation
    //Order parameter is the sum of distributions plus any source/correction terms

    if constexpr(std::tuple_size<typename traits::Forces>::value != 0){

        return CollisionBase<typename traits::Stencil>::computeZerothMoment(distribution)
        + std::apply([k](auto&... tests) {

            return (tests.computeDensitySource(k) + ...);

        }, ModelBase<traits>::mt_Forces);

    }
    else return CollisionBase<typename traits::Stencil>::computeZerothMoment(distribution);

}