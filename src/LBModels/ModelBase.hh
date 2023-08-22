#pragma once
#include <utility>
#include <memory>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include "../Forces/ForceBase.hh"
#include "../BoundaryModels/BoundaryBase.hh"
#include "../AddOns/AddOns.hh"
#include "../Collide.hh"
#include "../Service.hh"

template<class TLattice=void, int t_NumberOfComponents=1>
struct DefaultTrait : BaseTrait<DefaultTrait<TLattice,t_NumberOfComponents>> {
    
    using Stencil = std::conditional_t<TLattice::NDIM == 1, D1Q3, std::conditional_t<TLattice::NDIM == 2, D2Q9, D3Q19>>; //Here, D refers to the number of cartesian dimensions

    using Boundaries = std::tuple<>;

    using PreProcessors = std::tuple<>;

    using PostProcessors = std::tuple<>;

    using Forces = std::tuple<>;

    template<class TStencil>
    using CollisionModel = SRT<TStencil>;

    using Lattice = TLattice;

    static constexpr int NumberOfComponents = t_NumberOfComponents;

};

template<class TLattice, class TTraits = DefaultTrait<TLattice>>
class ModelBase { //Inherit from base class to avoid repetition of common
                                                      //calculations
    static_assert(std::is_base_of<StencilBase, typename TTraits::Stencil>(), "ERROR: invalid TStencil specified in TTraits class.");
    static_assert(TTraits::Stencil::D == TLattice::NDIM, "ERROR: The chosen TStencil must match the number of TLattice dimensions in the TLattice properties.");
    static_assert(CheckBaseTemplate<ForceBase, typename TTraits::Forces>::value, "ERROR: At least one TForce chosen is not a TForce class. The class must inherit from ForceBase.");
    static_assert(CheckBase<AddOnBase, typename TTraits::PreProcessors>::value, "ERROR: At least one TPreProcessor chosen is not an addon class. The class must inherit from AddOnBase.");
    static_assert(CheckBase<AddOnBase, typename TTraits::PostProcessors>::value, "ERROR: At least one TPostProcessor chosen is not an addon class.  The class must inherit from AddOnBase.");
    static_assert(CheckBase<BoundaryBase, typename TTraits::Boundaries>::value, "ERROR: At least one boundary condition chosen is not a boundary class. The class must inherit from BoundaryBase.");
    public:

        ModelBase()
            : mData(),
              mDistribution(mData.getDistributionObject())
        {
            TLattice latticeInit; // Initialise the TLattice and parallelisation
        }

        ModelBase(ModelBase<TLattice,TTraits>& other)
            : mData(other.mData),
              mDistribution(other.mDistribution)
        {}

        inline virtual void precompute(); //Perform any necessary computations before collision

        template<class TTupleType>
        inline void computeTuple(TTupleType& tup,int k); //Perform any necessary computations before collision

        template<class TTupleType>
        inline void precomputeTuple(TTupleType& tup,int k); //Perform any necessary computations before collision

        inline virtual void postprocess(); //Perform any necessary computations before collision

        template<class TTupleType>
        inline void postprocessTuple(TTupleType& tup,int k);

        template<class TTupleType>
        inline void communicateTuple(TTupleType& tup);

        template<class TTupleType>
        inline void communicatePrecomputeTuple(TTupleType& tup);

        template<class TTupleType>
        inline void communicatePostprocessTuple(TTupleType& tup);

        template<class TTupleType>
        inline void communicatePrecomputeBoundaries(TTupleType& tup);

        template<class TTupleType>
        inline void communicatePostprocessBoundaries(TTupleType& tup);

        inline virtual void collide() = 0; //Collision step

        inline virtual void boundaries(); //Boundary calculation

        inline virtual void initialise() = 0; //Initialisation step

        inline virtual void computeMomenta() = 0; //Momenta (density, velocity) calculation

        inline double computeForces(int xyz, int k) const; //Calculate other forces in direction xyz

        inline const std::vector<double>& getDistribution() const; //Return vector of distribution

        template<class TTupleType>
        inline auto getForceCalculator(TTupleType& TForceTuple, int k);

        template<class TForce>
        typename TForce::Method getMethod(TForce& f);

        template<class TForce, typename TForceTuple>
        void precomputeForces(TForce& f,TForceTuple& forcemethods,int k);

        template<class TPreProcessor, int inst = 0>
        inline TPreProcessor& getPreProcessor() {

            auto preprocessors = get_type<TPreProcessor>(mt_PreProcessors);

            return std::get<inst>(preprocessors);

        }

        template<class TPostProcessor, int inst = 0>
        inline TPostProcessor& getPostProcessor() {

            auto postprocessors = get_type<TPostProcessor>(mt_PostProcessors);

            return std::get<inst>(postprocessors);

        }

        template<class TForce, int inst = 0>
        inline TForce& getForce() {

            auto forces = get_type<TForce>(mt_Forces);

            return std::get<inst>(forces);

        }

        template<class boundary, int inst = 0>
        inline boundary& getBoundary() {

            auto boundaries = get_type<boundary>(mt_Boundaries);

            return std::get<inst>(boundaries);

        }

        template<class maptype, class prefactortuple, class forcetype>
        inline void setForceSums(prefactortuple& prefactors,forcetype& f, int idx, int k) {
            
            std::get<typename maptype::template get<typename forcetype::Prefactor>>(prefactors).val[idx] += f.template compute<TTraits>(idx, k);
            //std::cout<<forcesums[typeid(GuoPrefactor)][idx]<<std::endl;
        }

        inline void collisionQ(const double* equilibriums, const double* olddistributions, const double& inversetau, int k) {
            
            if constexpr(std::tuple_size<typename TTraits::Forces>::value != 0){

                if constexpr (mDistribution.SaveEquilibrium) {
                    mDistribution.saveEquilibriums(equilibriums,k);
                }
                
                auto forcemethods = getForceCalculator(mt_Forces,k);
                
                //std::unordered_map<std::type_index, std::array<double,TTraits::Stencil::Q>> forcesums;
                //std::remove_reference<decltype(forces)>::type::Method::Prefactor
                auto tempMap = std::apply([this](auto&... forces){//See Algorithm.hh for explanation of std::apply

                    ct_map_types<kv<typename std::remove_reference<decltype(forces)>::type::Method::Prefactor,std::array<double,TTraits::Stencil::Q>>...> tempmap;

                    return tempmap;

                }, mt_Forces);

                using ForcingMap = decltype(tempMap);

                auto tempTuple = std::apply([this](auto&... forces){//See Algorithm.hh for explanation of std::apply

                    
                    constexpr std::tuple<typename ForcingMap::template get<typename std::remove_reference<decltype(forces)>::type::Method::Prefactor>...> temptup;
                    constexpr auto temptup2 = make_tuple_unique(temptup);
                    return temptup2;

                }, mt_Forces);

                for (int idx = 0; idx < TTraits::Stencil::Q; idx++) {
                    std::apply([this,idx, k, &tempTuple](auto&... forces) mutable{
                        (this->setForceSums<ForcingMap>(tempTuple,forces, idx, k) , ...);
                        
                                                                    }, *forcemethods);
                    //std::cout<<forcesums[typeid(GuoPrefactor)][idx]<<std::endl;
                }
            
                //double i = forcesums[typeid(decltype(Guo::Prefactor))][0];
                auto tempforceprefactors=std::apply([](auto&... methods){//See Algorithm.hh for explanation of std::apply

                    return std::make_unique<decltype(make_tuple_unique(std::make_tuple(getForcePrefactor(methods))...))>();

                }, *forcemethods);

                for (int idx = 0; idx <TTraits::Stencil::Q; idx++) { //loop over discrete velocity directions
                //Set distribution at location "mDistribution.streamIndex" equal to the value returned by
                //"computeCollisionQ"
                
                    double collision=TTraits::template CollisionModel<typename TTraits::Stencil>::template collide<typename TTraits::Lattice>(olddistributions,equilibriums,inversetau,idx) 
                                     + std::apply([&inversetau,&tempTuple,idx,this](auto&... prefactors) mutable {
                                                        return (TTraits::template CollisionModel<typename TTraits::Stencil>::template forcing<typename TTraits::Lattice,decltype(prefactors)>(this->mt_Forces,&(std::get<typename ForcingMap::template get<typename remove_const_and_reference<decltype(prefactors)>::type>>(tempTuple).val[0]),inversetau,idx) + ...);
                                                                                                                }, *tempforceprefactors);
                    
                    mDistribution.getDistributionPointer(mDistribution.streamIndex(k, idx))[idx] = collision;

                }
            }
            else{
                for (int idx = 0; idx <TTraits::Stencil::Q; idx++) { //loop over discrete velocity directions
                //Set distribution at location "mDistribution.streamIndex" equal to the value returned by
                //"computeCollisionQ"
                
                    double collision=TTraits::template CollisionModel<typename TTraits::Stencil>::template collide<typename TTraits::Lattice>(olddistributions,equilibriums,inversetau,idx);

                    mDistribution.getDistributionPointer(mDistribution.streamIndex(k, idx))[idx] = collision;

                }
            }

        }

        template<class TForceTypes>
        void updateForces(double& TForce, TForceTypes& forcetypes, int k, int idx) {
            if constexpr(std::tuple_size<TForceTypes>::value != 0){

                TForce=std::apply([idx, k](auto&... forcetype){
                            return (forcetype.template compute<TTraits>(idx, k) + ...);
                                                                    }, forcetypes);

            }
        }

        inline double computeDensity(const double* distribution, const int k); //Calculate density

        inline double computeVelocity(const double* distribution, const double& density,
                                const int xyz, const int k); //Calculate velocity

        typename std::remove_reference<TLattice>::type::template DataType<typename TTraits::Stencil> mData; //MOVE THIS TO BASE
        typename std::remove_reference<TLattice>::type::template DataType<typename TTraits::Stencil>::DistributionData& mDistribution = mData.getDistributionObject();
            //Distributions

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions

        typename TTraits:: PreProcessors mt_PreProcessors; //MOVE THIS TO BASE
        typename TTraits:: PostProcessors mt_PostProcessors; //MOVE THIS TO BASE
        typename TTraits:: Forces mt_Forces; //MOVE THIS TO BASE
        typename TTraits:: Boundaries mt_Boundaries; //MOVE THIS TO BASE
        Geometry<TLattice> mGeometry; //MOVE THIS TO BASE

        
        std::vector<double>& distribution = mDistribution.getDistribution(); //Reference to vector of distributions
        
};

template<class TLattice, class TTraits>
inline double ModelBase<TLattice, TTraits>::computeDensity(const double* distribution, const int k) { //Density<> calculation
    //Density<> is the sum of distributions plus any source/correction terms

    if constexpr(std::tuple_size<typename TTraits::Forces>::value != 0) {

        return (CollisionBase<TLattice,typename TTraits::Stencil>::computeZerothMoment(distribution) + std::apply([k](auto&... forces) {

                return (forces.template computeDensitySource<TTraits>(k) + ...);

            }, mt_Forces)) * std::apply([k](auto&... forces) {

                return (forces.template computeDensitySourceMultiplicative<TTraits>(k) * ...);

            }, mt_Forces);

    }
    else return CollisionBase<TLattice,typename TTraits::Stencil>::computeZerothMoment(distribution);

}

template<class TLattice, class TTraits>
inline double ModelBase<TLattice, TTraits>::computeVelocity(const double* distribution, const double& density,
                                             const int xyz, const int k) { //Velocity calculation in direction xyz
    //Velocity in direction xyz is sum of distribution times the xyz component of the discrete velocity vector
    //in each direction plus any source/correction terms

    if constexpr(std::tuple_size<typename TTraits::Forces>::value != 0) {

        return (1./(Density<>::get<TLattice>(k)))*CollisionBase<TLattice,typename TTraits::Stencil>::computeFirstMoment(distribution, xyz) + (1./(Density<>::get<TLattice>(k)))*std::apply([xyz, k](auto&&... forces) mutable {

                return (forces.template computeVelocitySource<TTraits>(xyz, k) + ...);
                
            }, mt_Forces);

    }
    else return (1./(Density<>::get<TLattice>(k)))*CollisionBase<TLattice,typename TTraits::Stencil>::computeFirstMoment(distribution, xyz);

}

template<class TLattice, class TTraits>
inline const std::vector<double>& ModelBase<TLattice,TTraits>::getDistribution() const {

    return distribution; //Return reference to distribution vector

}

template<class TLattice, class TTraits>
template<class TForce, typename TForceTuple>
void ModelBase<TLattice,TTraits>::precomputeForces(TForce& f,TForceTuple& forcemethods,int k){
    std::get<decltype(getMethod(f))>(forcemethods).template precompute<TTraits>(f,k);
}

template<class TLattice, class TTraits>
template<class TTupleType>
inline auto ModelBase<TLattice,TTraits>::getForceCalculator(TTupleType& TForceTuple, int k){

    if constexpr(std::tuple_size<TTupleType>::value != 0){
        
        auto tempforce=std::apply([k](auto&... forces){//See Algorithm.hh for explanation of std::apply

            return std::make_unique<decltype(make_tuple_unique(std::make_tuple(getMethod(forces))...))>();

        }, TForceTuple);

        tempforce=std::apply([this,tempforce=std::move(tempforce),k](auto&... forces) mutable {//See Algorithm.hh for explanation of std::apply
            
            (this->precomputeForces(forces,*tempforce,k),...);
            return std::move(tempforce);
        
        }, TForceTuple);

        return tempforce;

    }
    else return std::make_unique<std::tuple<>>();

}

template<class TLattice, class TTraits>
inline void ModelBase<TLattice,TTraits>::precompute() {

    #pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k <TLattice::N - TLattice::HaloSize; k++) { //loop over k

        precomputeTuple(mt_Boundaries,k);

        computeTuple(mt_PreProcessors,k);

        precomputeTuple(mt_Forces,k);
        
    }
  
    communicateTuple(mt_PreProcessors);
    communicatePrecomputeTuple(mt_Forces);
    communicatePrecomputeBoundaries(mt_Boundaries);

    #pragma omp master
    {
    mDistribution.getDistribution().swap(mDistribution.getDistributionOld()); //swap old and new distributions
                                                                                //before collision
    }
    
}

template<class TLattice, class TTraits>
template<class TTupleType>
inline void ModelBase<TLattice,TTraits>::computeTuple(TTupleType& tup, int k) {

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([k](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template compute<TTraits>(k),...);

        }, tup);

    }
    
}

template<class TLattice, class TTraits>
template<class TTupleType>
inline void ModelBase<TLattice,TTraits>::precomputeTuple(TTupleType& tup, int k) {

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([k](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template precompute<TTraits>(k),...);

        }, tup);

    }
    
}

template<class TLattice, class TTraits>
template<class TTupleType>
inline void ModelBase<TLattice,TTraits>::postprocessTuple(TTupleType& tup, int k) {

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([k](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template postprocess<TTraits>(k),...);

        }, tup);

    }
    
}

template<class TLattice, class TTraits>
template<class TTupleType>
inline void ModelBase<TLattice,TTraits>::communicateTuple(TTupleType& tup) {

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template communicate<TTraits>(),...);

        }, tup);

    }
    
}

template<class TLattice, class TTraits>
template<class TTupleType>
inline void ModelBase<TLattice,TTraits>::communicatePrecomputeTuple(TTupleType& tup) {

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template communicatePrecompute<TTraits>(),...);

        }, tup);

    }
    
}

template<class TLattice, class TTraits>
template<class TTupleType>
inline void ModelBase<TLattice,TTraits>::communicatePostprocessTuple(TTupleType& tup) {

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template communicatePostProcess<TTraits>(),...);

        }, tup);

    }
    
}

template<class TLattice, class TTraits>
template<class TTupleType>
inline void ModelBase<TLattice,TTraits>::communicatePrecomputeBoundaries(TTupleType& tup) {

    communicatePrecomputeTuple(mt_Boundaries);

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([this](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template communicatePrecompute<TTraits>(this -> mDistribution),...);

        }, tup);

    }
    
}

template<class TLattice, class TTraits>
template<class TTupleType>
inline void ModelBase<TLattice,TTraits>::communicatePostprocessBoundaries(TTupleType& tup) {

    communicatePostprocessTuple(mt_Boundaries);

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([this](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template communicatePostProcess<TTraits>(this -> mDistribution),...);

        }, tup);

    }
    
}

template<class TLattice, class TTraits>
inline void ModelBase<TLattice,TTraits>::postprocess() {

    #pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k <TLattice::N - TLattice::HaloSize; k++) { //loop over k

        postprocessTuple(mt_Boundaries,k);

        computeTuple(mt_PostProcessors,k);

        postprocessTuple(mt_Forces,k);
        
    }
  
    communicateTuple(mt_PostProcessors);
    communicatePostprocessTuple(mt_Forces);
    communicatePrecomputeBoundaries(mt_Boundaries);
    
}

template<class TLattice, class TTraits>
inline void ModelBase<TLattice,TTraits>::boundaries() {

    #pragma omp for schedule(guided)
    for (int k = 0; k <TLattice::N; k++) { //loop over k

        if constexpr(std::tuple_size<typename TTraits::Boundaries>::value != 0) { //Check if there are any boundary
                                                                              //models
            std::apply([this, k](auto&... boundaries) {
                        
                                (boundaries.template compute<TTraits>(this -> mDistribution, k) , ...); 

                    }, mt_Boundaries);
            
        }

    }
    
}
