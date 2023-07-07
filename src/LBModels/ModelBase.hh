#pragma once
#include <utility>
#include <memory>
#include "../Forces/ForceBase.hh"
#include "../BoundaryModels/BoundaryBase.hh"
#include "../AddOns/AddOns.hh"
#include "../Collide.hh"
#include "../Service.hh"

template<class t_Lattice=void, int t_NumberOfComponents=1>
struct DefaultTrait : BaseTrait<DefaultTrait<t_Lattice,t_NumberOfComponents>> {
    
    using Stencil = std::conditional_t<t_Lattice::NDIM == 1, D1Q3, std::conditional_t<t_Lattice::NDIM == 2, D2Q9, D3Q19>>; //Here, D refers to the number of cartesian dimensions

    using Boundaries = std::tuple<>;

    using PreProcessors = std::tuple<>;

    using PostProcessors = std::tuple<>;

    using Forces = std::tuple<>;

    template<class stencil>
    using CollisionModel = SRT<stencil>;

    using Lattice = t_Lattice;

    static constexpr int NumberOfComponents = t_NumberOfComponents;

};

template<class lattice, class traits = DefaultTrait<lattice>>
class ModelBase { //Inherit from base class to avoid repetition of common
                                                      //calculations
    static_assert(std::is_base_of<StencilBase, typename traits::Stencil>(), "ERROR: invalid stencil specified in traits class.");
    static_assert(traits::Stencil::D == lattice::NDIM, "ERROR: The chosen stencil must match the number of lattice dimensions in the lattice properties.");
    static_assert(CheckBaseTemplate<ForceBase, typename traits::Forces>::value, "ERROR: At least one boundary condition chosen is not a boundary class.");
    static_assert(CheckBase<AddOnBase, typename traits::PreProcessors>::value, "ERROR: At least one boundary condition chosen is not a boundary class.");
    static_assert(CheckBase<AddOnBase, typename traits::PostProcessors>::value, "ERROR: At least one boundary condition chosen is not a boundary class.");
    static_assert(CheckBase<BoundaryBase, typename traits::Boundaries>::value, "ERROR: At least one boundary chosen is not a boundary class.");
    public:

        ModelBase()
            : m_Data(),
              m_Distribution(m_Data.getDistributionObject())
        {
            lattice latticeInit; // Initialise the lattice and parallelisation
        }

        ModelBase(ModelBase<lattice,traits>& other)
            : m_Data(other.m_Data),
              m_Distribution(other.m_Distribution)
        {}

        inline virtual void precompute(); //Perform any necessary computations before collision

        template<class tupletype>
        inline void computeTuple(tupletype& tup,int k); //Perform any necessary computations before collision

        template<class tupletype>
        inline void precomputeTuple(tupletype& tup,int k); //Perform any necessary computations before collision

        inline virtual void postprocess(); //Perform any necessary computations before collision

        template<class tupletype>
        inline void postprocessTuple(tupletype& tup,int k);

        template<class tupletype>
        inline void communicateTuple(tupletype& tup);

        template<class tupletype>
        inline void communicatePrecomputeTuple(tupletype& tup);

        template<class tupletype>
        inline void communicatePostprocessTuple(tupletype& tup);

        inline virtual void collide() = 0; //Collision step

        inline virtual void boundaries(); //Boundary calculation

        inline virtual void initialise() = 0; //Initialisation step

        inline virtual void computeMomenta() = 0; //Momenta (density, velocity) calculation

        inline double computeForces(int xyz, int k) const; //Calculate other forces in direction xyz

        inline const std::vector<double>& getDistribution() const; //Return vector of distribution

        template<class tupletype>
        inline auto getForceCalculator(tupletype& forcetuple, int k);

        template<class force>
        typename force::Method getMethod(force& f);

        template<class force, class forcetuple>
        void precomputeForces(force& f,forcetuple& forcemethods,int k);

        template<class preprocessor, int inst = 0>
        inline preprocessor& getPreProcessor() {

            auto preprocessors = get_type<preprocessor>(mt_PreProcessors);

            return std::get<inst>(preprocessors);

        }

        template<class postprocessor, int inst = 0>
        inline postprocessor& getPostProcessor() {

            auto postprocessors = get_type<postprocessor>(mt_PostProcessors);

            return std::get<inst>(postprocessors);

        }

        template<class force, int inst = 0>
        inline force& getForce() {

            auto forces = get_type<force>(mt_Forces);

            return std::get<inst>(forces);

        }

        template<class boundary, int inst = 0>
        inline boundary& getBoundary() {

            auto boundaries = get_type<boundary>(mt_Boundaries);

            return std::get<inst>(boundaries);

        }

        inline void collisionQ(const double* forces, const double* equilibriums, const double* olddistributions, const double& inversetau, int k) {

            for (int idx = 0; idx <traits::Stencil::Q; idx++) { //loop over discrete velocity directions
                //Set distribution at location "m_Distribution.streamIndex" equal to the value returned by
                //"computeCollisionQ"
                
                double collision=traits::template CollisionModel<typename traits::Stencil>::template collide<typename traits::Lattice>(olddistributions,equilibriums,inversetau,idx);

                if constexpr(std::tuple_size<typename traits::Forces>::value != 0){
                    collision += traits::template CollisionModel<typename traits::Stencil>::template forcing<typename traits::Lattice>(forces,inversetau,idx);
                    
                }

                m_Distribution.getDistributionPointer(m_Distribution.streamIndex(k, idx))[idx] = collision;

            }

        }

        template<class T_forcetypes>
        void updateForces(double& force, T_forcetypes& forcetypes, int k, int idx) {
            if constexpr(std::tuple_size<T_forcetypes>::value != 0){

                force=std::apply([idx, k](auto&... forcetype){
                            return (forcetype.template compute<traits>(idx, k) + ...);
                                                                    }, forcetypes);

            }
        }

        inline double computeDensity(const double* distribution, const int k); //Calculate density

        inline double computeVelocity(const double* distribution, const double& density,
                                const int xyz, const int k); //Calculate velocity

        typename std::remove_reference<lattice>::type::template DataType<typename traits::Stencil> m_Data; //MOVE THIS TO BASE
        typename std::remove_reference<lattice>::type::template DataType<typename traits::Stencil>::DistributionData& m_Distribution = m_Data.getDistributionObject();
            //Distributions

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions

        typename traits:: PreProcessors mt_PreProcessors; //MOVE THIS TO BASE
        typename traits:: PostProcessors mt_PostProcessors; //MOVE THIS TO BASE
        typename traits:: Forces mt_Forces; //MOVE THIS TO BASE
        typename traits:: Boundaries mt_Boundaries; //MOVE THIS TO BASE
        Geometry<lattice> m_Geometry; //MOVE THIS TO BASE

        
        std::vector<double>& distribution = m_Distribution.getDistribution(); //Reference to vector of distributions
        
};

template<class T_lattice, class T_traits>
inline double ModelBase<T_lattice, T_traits>::computeDensity(const double* distribution, const int k) { //Density<> calculation
    //Density<> is the sum of distributions plus any source/correction terms

    if constexpr(std::tuple_size<typename T_traits::Forces>::value != 0) {

        return CollisionBase<T_lattice,typename T_traits::Stencil>::computeZerothMoment(distribution) + std::apply([k](auto&... forces) {

                return (forces.template computeDensitySource<T_traits>(k) + ...);

            }, mt_Forces);

    }
    else return CollisionBase<T_lattice,typename T_traits::Stencil>::computeZerothMoment(distribution);

}

template<int i>
void test(){}

template<class T_lattice, class T_traits>
inline double ModelBase<T_lattice, T_traits>::computeVelocity(const double* distribution, const double& density,
                                             const int xyz, const int k) { //Velocity calculation in direction xyz
    //Velocity in direction xyz is sum of distribution times the xyz component of the discrete velocity vector
    //in each direction plus any source/correction terms

    if constexpr(std::tuple_size<typename T_traits::Forces>::value != 0) {
        std::get<BodyForce<>>(mt_Forces);
        //test<decltype(std::get<BodyForce>(mt_Forces))>();
        return (1./(Density<>::get<T_lattice>(k)))*CollisionBase<T_lattice,typename T_traits::Stencil>::computeFirstMoment(distribution, xyz) + (1./(Density<>::get<T_lattice>(k)))*std::apply([xyz, k](auto&&... forces) mutable {

                //return (test<decltype(forces)>() + ...);
                return (forces.template computeVelocitySource<T_traits>(xyz, k) + ...);
                
            }, mt_Forces);

    }
    else return (1./(Density<>::get<T_lattice>(k)))*CollisionBase<T_lattice,typename T_traits::Stencil>::computeFirstMoment(distribution, xyz);

}

template<class lattice, class traits>
inline const std::vector<double>& ModelBase<lattice,traits>::getDistribution() const {

    return distribution; //Return reference to distribution vector

}

template<class lattice, class traits>
template<class force>
typename force::Method ModelBase<lattice,traits>::getMethod(force& f){
    return std::declval<typename force::Method>();
}

template<class lattice, class traits>
template<class force, class forcetuple>
void ModelBase<lattice,traits>::precomputeForces(force& f,forcetuple& forcemethods,int k){
    std::get<decltype(getMethod(f))>(forcemethods).template precompute<traits>(f,k);
}

template<class lattice, class traits>
template<class tupletype>
inline auto ModelBase<lattice,traits>::getForceCalculator(tupletype& forcetuple, int k){

    if constexpr(std::tuple_size<tupletype>::value != 0){
        
        auto tempforce=std::apply([k](auto&... forces){//See Algorithm.hh for explanation of std::apply

            return std::make_unique<decltype(make_tuple_unique(std::make_tuple(getMethod(forces))...))>();

        }, forcetuple);

        tempforce=std::apply([this,tempforce=std::move(tempforce),k](auto&... forces) mutable {//See Algorithm.hh for explanation of std::apply
            
            (this->precomputeForces(forces,*tempforce,k),...);
            return std::move(tempforce);
        
        }, forcetuple);

        return tempforce;

    }
    else return std::make_unique<std::tuple<>>();

}

template<class lattice, class traits>
inline void ModelBase<lattice,traits>::precompute() {

    #pragma omp for schedule(guided)
    for (int k = lattice::HaloSize; k <lattice::N - lattice::HaloSize; k++) { //loop over k

        precomputeTuple(mt_Boundaries,k);

        computeTuple(mt_PreProcessors,k);

        precomputeTuple(mt_Forces,k);
        
    }
  
    communicateTuple(mt_PreProcessors);
    communicatePrecomputeTuple(mt_Forces);

    #pragma omp master
    {
    m_Distribution.getDistribution().swap(m_Distribution.getDistributionOld()); //swap old and new distributions
                                                                                //before collision
    }
    
}

template<class lattice, class traits>
template<class tupletype>
inline void ModelBase<lattice,traits>::computeTuple(tupletype& tup, int k) {

    if constexpr(std::tuple_size<tupletype>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([k](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template compute<traits>(k),...);

        }, tup);

    }
    
}

template<class lattice, class traits>
template<class tupletype>
inline void ModelBase<lattice,traits>::precomputeTuple(tupletype& tup, int k) {

    if constexpr(std::tuple_size<tupletype>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([k](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template precompute<traits>(k),...);

        }, tup);

    }
    
}

template<class lattice, class traits>
template<class tupletype>
inline void ModelBase<lattice,traits>::postprocessTuple(tupletype& tup, int k) {

    if constexpr(std::tuple_size<tupletype>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([k](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template postprocess<traits>(k),...);

        }, tup);

    }
    
}

template<class lattice, class traits>
template<class tupletype>
inline void ModelBase<lattice,traits>::communicateTuple(tupletype& tup) {

    if constexpr(std::tuple_size<tupletype>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template communicate<traits>(),...);

        }, tup);

    }
    
}

template<class lattice, class traits>
template<class tupletype>
inline void ModelBase<lattice,traits>::communicatePrecomputeTuple(tupletype& tup) {

    if constexpr(std::tuple_size<tupletype>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template communicatePrecompute<traits>(),...);

        }, tup);

    }
    
}

template<class lattice, class traits>
template<class tupletype>
inline void ModelBase<lattice,traits>::communicatePostprocessTuple(tupletype& tup) {

    if constexpr(std::tuple_size<tupletype>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template communicatePostProcess<traits>(),...);

        }, tup);

    }
    
}

template<class lattice, class traits>
inline void ModelBase<lattice,traits>::postprocess() {

    #pragma omp for schedule(guided)
    for (int k = lattice::HaloSize; k <lattice::N - lattice::HaloSize; k++) { //loop over k

        postprocessTuple(mt_Boundaries,k);

        computeTuple(mt_PostProcessors,k);

        postprocessTuple(mt_Forces,k);
        
    }
  
    communicateTuple(mt_PostProcessors);
    communicatePostprocessTuple(mt_Forces);
    
}

template<class lattice, class traits>
inline void ModelBase<lattice,traits>::boundaries() {

    #pragma omp for schedule(guided)
    for (int k = 0; k <lattice::N; k++) { //loop over k

        if constexpr(std::tuple_size<typename traits::Boundaries>::value != 0) { //Check if there are any boundary
                                                                              //models
            for (int idx = 0; idx <traits::Stencil::Q; idx++) {

                if(m_Geometry.isSolid(k) && !m_Geometry.isSolid(m_Distribution.streamIndex(k, idx))) {
                    std::apply([this, k, idx](auto&... boundaries) {
                        // Make this a sub function for readability
                                (boundaries.template compute<traits>(this -> m_Distribution, k, idx) , ...);

                    }, mt_Boundaries);
                }

            }
            
        }

    }
    
}
