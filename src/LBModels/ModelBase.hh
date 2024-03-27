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
#include "../Data.hh"

template<class TLattice=void, int t_NumberOfComponents=1>
struct DefaultTrait : BaseTrait<DefaultTrait<TLattice,t_NumberOfComponents>> {
    
    using Stencil = std::conditional_t<TLattice::NDIM == 1, D1Q3, std::conditional_t<TLattice::NDIM == 2, D2Q9, D3Q19>>; //Here, D refers to the number of cartesian dimensions

    using Boundaries = std::tuple<std::tuple<>>;

    using Processors = std::tuple<>;

    using Forces = std::tuple<>;

    template<class TStencil>
    using CollisionModel = SRT<TStencil>;

    using Lattice = TLattice;

    template<class Tlattice, class TStencil>
    using DataType = DataOldNew<TLattice,TStencil>;

    static constexpr int NumberOfComponents = t_NumberOfComponents;

};


class Model {}; // Used to store the model without template parameters


template<class TLattice, class TTraits = DefaultTrait<TLattice>>
class ModelBase : public Model {
    static_assert(std::is_base_of<StencilBase, typename TTraits::Stencil>(), "ERROR: invalid TStencil specified in TTraits class.");
    static_assert(TTraits::Stencil::D == TLattice::NDIM, "ERROR: The chosen TStencil must match the number of TLattice dimensions in the TLattice properties.");
    static_assert(CheckBaseTemplate<ForceBase, typename TTraits::Forces>::value, "ERROR: At least one TForce chosen is not a TForce class. The class must inherit from ForceBase.");
    //static_assert(CheckBase<AddOnBase, typename TTraits::Processors>::value, "ERROR: At least one TPreProcessor chosen is not an addon class. The class must inherit from AddOnBase.");
    //static_assert(CheckBase<AddOnBase, typename TTraits::PostProcessors>::value, "ERROR: At least one TPostProcessor chosen is not an addon class.  The class must inherit from AddOnBase.");
    //static_assert(CheckBase<BoundaryBase, typename TTraits::Boundaries>::value, "ERROR: At least one boundary condition chosen is not a boundary class. The class must inherit from BoundaryBase.");
    public:

        ModelBase()
            : latticeInit(),
              mData(),
              mDistribution(mData.getDistributionObject())
        {
             
             // Initialise the TLattice and parallelisation
        }

        ModelBase(ModelBase<TLattice,TTraits>& other)
            : latticeInit(),
              mData(other.mData),
              mDistribution(other.mDistribution)
        {}
        /*
        inline virtual void precompute(); //Perform any necessary computations before collision
        */
        template<class TTupleType>
        inline void computeTuple(TTupleType& tup,int k); //Perform any necessary computations before collision
        /*
        template<class TTupleType>
        inline void precomputeTuple(TTupleType& tup,int k); //Perform any necessary computations before collision
        */
        inline virtual void computeProcessors(); //Perform any necessary computations before collision

        template<class TTupleType>
        inline void processorTuple(TTupleType& tup,int k);

        template<class TTupleType>
        inline void communicateTuple(TTupleType& tup);

        template<class TTupleType>
        inline void communicateProcessorBoundaries(TTupleType& tup);
        /*
        template<class TTupleType>
        inline void communicatePrecomputeTuple(TTupleType& tup);

        template<class TTupleType>
        inline void communicatePostprocessTuple(TTupleType& tup);

        template<class TTupleType>
        inline void communicatePrecomputeBoundaries(TTupleType& tup);

        template<class TTupleType>
        inline void communicatePostprocessBoundaries(TTupleType& tup);
        */
        inline virtual void collide() = 0; //Collision step

        inline virtual void stream(); //Collision step

        inline virtual void boundaries(); //Boundary calculation

        inline bool isNan() {
            bool isnan = false;
            #pragma omp for schedule(guided)
            for (auto i : distribution)
            {
                if (std::isnan(i)||std::isinf(i)) isnan = true;
            }
            return isnan;
        }

        template<class TBoundaryType>
        inline void runBoundaries(TBoundaryType& boundary); //Boundary calculation

        template<class TProcessorType>
        inline void runProcessors(TProcessorType& processor); //Boundary calculation

        inline virtual void initialise() = 0; //Initialisation step

        inline virtual void computeMomenta() = 0; //Momenta (density, velocity) calculation

        inline virtual double computeEquilibrium(int k, int idx) = 0; //Calculate equilibrium in direction idx

        inline virtual double computePressure(int k); //Pressure calculation

        inline double computeForces(int xyz, int k) const; //Calculate other forces in direction xyz

        inline const std::vector<double>& getDistribution() const; //Return vector of distribution

        inline void initialiseBoundaries();

        template<class TTupleType>
        inline auto getForceCalculator(TTupleType& TForceTuple, int k);

        template<class TForce>
        typename TForce::Method getMethod(TForce& f);

        template<class TForce, typename TForceTuple>
        void precomputeForces(TForce& f,TForceTuple& forcemethods,int k);

        template<class TProcessor, int tuplenum, int inst=0>
        inline TProcessor& getProcessor() {
            static_assert(has_type<TProcessor,std::tuple_element_t<tuplenum,typename TTraits::Processors>>::value, "Desired processor is not included in the model.");

            auto processors = get_type<TProcessor>(std::get<tuplenum>(mt_Processors));

            return std::get<inst>(processors);

        }

        template<class TProcessor>
        inline TProcessor& getProcessor() {

            using indices = typename tuple_tuple_index_of<typename TTraits::Processors,TProcessor>::idx;

            static constexpr int idx1 = std::tuple_element<0,indices>::type::value;
            static constexpr int idx2 = std::tuple_element<1,indices>::type::value;

            return std::get<idx2>(std::get<idx1>(mt_Processors));

        }

        template<class TForce, int inst = 0>
        inline TForce& getForce() {
            static_assert(has_type<TForce,typename TTraits::Forces>::value, "Desired force is not included in the model.");

            auto forces = get_type<TForce>(mt_Forces);

            return std::get<inst>(forces);

        }

        template<class TBoundary, int tuplenum, int inst=0>
        inline TBoundary& getBoundary() {
            static_assert(has_type<TBoundary,std::tuple_element_t<tuplenum,typename TTraits::Boundaries>>::value, "Desired boundary is not included in the model.");

            auto boundaries = get_type<TBoundary>(std::get<tuplenum>(mt_Boundaries));

            return std::get<inst>(boundaries);

        }

        template<class TBoundary>
        inline TBoundary& getBoundary() {

            using indices = typename tuple_tuple_index_of<typename TTraits::Boundaries,TBoundary>::idx;
            
            static constexpr int idx1 = std::tuple_element<0,indices>::type::value;
            static constexpr int idx2 = std::tuple_element<1,indices>::type::value;

            return std::get<idx2>(std::get<idx1>(mt_Boundaries));

        }

        template<class maptype, class prefactortuple, class forcetype>
        inline void setForceSums(prefactortuple& prefactors,forcetype& f, int idx, int k) {
            
            std::get<typename maptype::template get<typename forcetype::Prefactor>>(prefactors).val[idx] += f.template compute<TTraits>(idx, k);
            //std::cout<<forcesums[typeid(GuoPrefactor)][idx]<<std::endl;
        }

        inline void collisionQ(const double* equilibriums, const double* olddistributions, const double& inversetau, int k) {
            
            if constexpr(std::tuple_size<typename TTraits::Forces>::value != 0){

                if (mDistribution.SaveEquilibrium) {
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
                #pragma omp simd
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
                #pragma omp simd
                for (int idx = 0; idx <TTraits::Stencil::Q; idx++) { //loop over discrete velocity directions
                //Set distribution at location "mDistribution.streamIndex" equal to the value returned by
                //"computeCollisionQ"
                
                    double collision=TTraits::template CollisionModel<typename TTraits::Stencil>::template collide<typename TTraits::Lattice>(olddistributions,equilibriums,inversetau,idx) 
                                     + std::apply([&inversetau,&tempTuple,idx,this](auto&... prefactors) mutable {
                                                        return (TTraits::template CollisionModel<typename TTraits::Stencil>::template forcing<typename TTraits::Lattice,decltype(prefactors)>(this->mt_Forces,&(std::get<typename ForcingMap::template get<typename remove_const_and_reference<decltype(prefactors)>::type>>(tempTuple).val[0]),inversetau,idx) + ...);
                                                                                                                }, *tempforceprefactors);
                    
                    mDistribution.getPostCollisionDistribution(k,idx) = collision;

                }
            }
            else{
                
                for (int idx = 0; idx <TTraits::Stencil::Q; idx++) { //loop over discrete velocity directions
                //Set distribution at location "mDistribution.streamIndex" equal to the value returned by
                //"computeCollisionQ"
                
                    double collision=TTraits::template CollisionModel<typename TTraits::Stencil>::template collide<typename TTraits::Lattice>(olddistributions,equilibriums,inversetau,idx);
                    
                    mDistribution.getPostCollisionDistribution(k,idx) = collision;

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

        static inline double computeVelocity(const double* distribution, typename TTraits::Forces& forcetuple, const double& density,
                                const int xyz, const int k); //Calculate velocity

        TLattice latticeInit;
        typename TTraits::template DataType<TLattice,typename TTraits::Stencil> mData;
        typename TTraits::template DataType<TLattice,typename TTraits::Stencil>::DistributionData& mDistribution = mData.getDistributionObject();
            //Distributions

        enum{ x = 0, y = 1, z = 2 }; //Indices corresponding to x, y, z directions

        typename TTraits:: Processors mt_Processors;
        typename TTraits:: Forces mt_Forces;
        typename TTraits:: Boundaries mt_Boundaries;
        Geometry<TLattice> mGeometry;

        inline void setCollideID(int id) {mCollideIDs[0]=id;};

        inline void setCollideID(const std::vector<int>& id) {mCollideIDs=id;};

        inline bool isCollisionNode(int k) {
            for (int i : mCollideIDs){
                if(Geometry<TLattice>::getBoundaryType(k) == i) return true;
            }
            return false;
        }
        std::vector<int> mCollideIDs = {0};
        
        std::vector<double>& distribution = mDistribution.getDistribution(); //Reference to vector of distributions
        
};


template<class TLattice, class TTraits>
inline void ModelBase<TLattice, TTraits>::initialiseBoundaries() {
    std::apply([this](auto&... boundarytuple){(std::apply([this](auto&... boundary){(boundary.initialise(this), ...);}, boundarytuple),...);}, mt_Boundaries);
}


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
inline double ModelBase<TLattice, TTraits>::computeVelocity(const double* distribution, typename TTraits::Forces& forcetuple, const double& density,
                                             const int xyz, const int k) { //Velocity calculation in direction xyz
    //Velocity in direction xyz is sum of distribution times the xyz component of the discrete velocity vector
    //in each direction plus any source/correction terms

    if constexpr(std::tuple_size<typename TTraits::Forces>::value != 0) {

        return (1./(Density<>::get<TLattice>(k)))*CollisionBase<TLattice,typename TTraits::Stencil>::computeFirstMoment(distribution, xyz) + (1./(Density<>::get<TLattice>(k)))*std::apply([xyz, k](auto&&... forces) mutable {

                return (forces.template computeVelocitySource<TTraits>(xyz, k) + ...);
                
            }, forcetuple);

    }
    else return (1./(Density<>::get<TLattice>(k)))*CollisionBase<TLattice,typename TTraits::Stencil>::computeFirstMoment(distribution, xyz);

}


template<class TLattice, class TTraits>
inline double ModelBase<TLattice, TTraits>::computePressure(int k) {
    double density = computeDensity(this->mDistribution.getDistributionPointer(k), k);
    return density * TTraits::Stencil::Cs2;
}


template<class TLattice, class TTraits>
inline const std::vector<double>& ModelBase<TLattice,TTraits>::getDistribution() const {

    return distribution; //Return reference to distribution vector

}

template<class TLattice, class TTraits>
template<class TTupleType>
inline auto ModelBase<TLattice,TTraits>::getForceCalculator(TTupleType& TForceTuple, int k){

    if constexpr(std::tuple_size<TTupleType>::value != 0){
        
        auto tempforce=std::apply([k,this](auto&... forces){//See Algorithm.hh for explanation of std::apply

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
inline void ModelBase<TLattice,TTraits>::stream() {

    TLattice::ResetParallelTracking();

    if constexpr (TTraits::template DataType<TLattice,typename TTraits::Stencil>::IsStreamingSeperate){

        #pragma omp for schedule(guided)
        for (int k = TLattice::HaloSize; k <TLattice::N - TLattice::HaloSize; k++) { //loop over k

            for (int idx = 0; idx <TTraits::Stencil::Q; idx++) {

                mDistribution.getPostStreamingDistribution(k,idx) = mDistribution.getPostCollisionDistribution(k)[idx]; 
            
            }

        }

    }

    this -> mData.communicateDistribution();
    std::apply([this](auto&... boundaryprocessor) {
            (communicateTuple(boundaryprocessor),...);
        }, mt_Boundaries);

    TLattice::ResetParallelTracking();
    
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
inline void ModelBase<TLattice,TTraits>::communicateProcessorBoundaries(TTupleType& tup) {

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template communicateProcessor<TTraits>(),...);

        }, tup);

    }
    
}
/*
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

    //communicatePrecomputeTuple(mt_Boundaries);

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([this](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template communicatePrecompute<TTraits>(),...);

        }, tup);

    }
    
}

template<class TLattice, class TTraits>
template<class TTupleType>
inline void ModelBase<TLattice,TTraits>::communicatePostprocessBoundaries(TTupleType& tup) {

    //communicatePostprocessTuple(mt_Boundaries);

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([this](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template communicatePostProcess<TTraits>(),...);

        }, tup);

    }
    
}
*/
template<class TLattice, class TTraits>
template<class TForce, typename TForceTuple>
void ModelBase<TLattice,TTraits>::precomputeForces(TForce& f,TForceTuple& forcemethods,int k){
    std::get<decltype(getMethod(f))>(forcemethods).template precompute<TTraits>(f,k);
}
/*
template<class TLattice, class TTraits>
inline void ModelBase<TLattice,TTraits>::postprocess() {

    TLattice::ResetParallelTracking();

    #pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k <TLattice::N - TLattice::HaloSize; k++) { //loop over k
        
        std::apply([this,k](auto&... boundaryprocessor) {
            (postprocessTuple(boundaryprocessor,k),...);
        }, mt_Boundaries);

        computeTuple(mt_PostProcessors,k);

        postprocessTuple(mt_Forces,k);
        
    }
  
    communicateTuple(mt_PostProcessors);
    communicatePostprocessTuple(mt_Forces);
    std::apply([this](auto&... boundaryprocessor) {
            (communicatePostprocessBoundaries(boundaryprocessor),...);
        }, mt_Boundaries);

    TLattice::ResetParallelTracking();
    
}
*/

template<class TLattice, class TTraits>
template<class TTupleType>
inline void ModelBase<TLattice,TTraits>::processorTuple(TTupleType& tup, int k) {

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([k](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template runProcessor<TTraits>(k),...);

        }, tup);

    }
    
}

template<class TLattice, class TTraits>
inline void ModelBase<TLattice,TTraits>::computeProcessors() {

    TLattice::ResetParallelTracking();

    std::apply([this](auto&... processor) {

        (runProcessors(processor), ...);
                
    }, mt_Processors);

    std::apply([this](auto&... boundaryprocessor) {
            (communicateProcessorBoundaries(boundaryprocessor),...);
        }, mt_Boundaries);

    #pragma omp for schedule(guided)
    for (int k = TLattice::HaloSize; k <TLattice::N - TLattice::HaloSize; k++) { //loop over k

        std::apply([this,k](auto&... boundaryprocessor) {
            (processorTuple(boundaryprocessor,k),...);
        }, mt_Boundaries);

        processorTuple(mt_Forces,k);
        
    }

    communicateTuple(mt_Forces);
    

    TLattice::ResetParallelTracking();

    #pragma omp master
    {
    mDistribution.getDistribution().swap(mDistribution.getDistributionOld()); //swap old and new distributions
                                                                                //before collision
    }
    
}

/*
template<class TLattice, class TTraits>
template<class TTupleType>
inline void ModelBase<TLattice,TTraits>::processorTuple(TTupleType& tup, int k) {

    if constexpr(std::tuple_size<TTupleType>::value != 0){ //Check if there is at least one element
                                                                        //in F

        std::apply([k](auto&... obj){//See Algorithm.hh for explanation of std::apply

            (obj.template runProcessor<TTraits>(k),...);

        }, tup);

    }
    
}*/
/*
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
*/
template<class TLattice, class TTraits>
inline void ModelBase<TLattice,TTraits>::boundaries() {

    TLattice::ResetParallelTracking();
    
    std::apply([this](auto&... boundaryprocessor) {

        (runBoundaries(boundaryprocessor), ...);
                
        }, mt_Boundaries);

    //this -> mData.communicateDistribution();
    
    TLattice::ResetParallelTracking();
    
}

template<class TLattice, class TTraits>
template<class TBoundaryType>
inline void ModelBase<TLattice,TTraits>::runBoundaries(TBoundaryType& boundaryprocessor) {

    if constexpr(std::tuple_size<typename TTraits::Boundaries>::value != 0) {

        #pragma omp for schedule(guided)
        for (int k = 0; k <TLattice::N; k++) { //loop over k
            //Check if there are any boundary
                                                                                //models
            std::apply([this, k](auto&... boundaries) {
                        
                                (boundaries.template compute<TTraits>(this -> mDistribution, k) , ...); 

                    }, boundaryprocessor);

        }

    }
}

template<class TLattice, class TTraits>
template<class TProcessorType>
inline void ModelBase<TLattice,TTraits>::runProcessors(TProcessorType& processortuple) {

    if constexpr(std::tuple_size<typename TTraits::Processors>::value != 0) {

        #pragma omp for schedule(guided)
        for (int k = 0; k <TLattice::N; k++) { //loop over k
            
            std::apply([this, k](auto&... processor) {
                        
                        (processor.template compute<TTraits>(k) , ...); 

                    }, processortuple);

        }
        communicateTuple(processortuple);
    }
}
