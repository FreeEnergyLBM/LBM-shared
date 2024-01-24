#pragma once
#include "../Geometry.hh"

class Model;


class BoundaryBase{
    public:

        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits>
        inline void communicateProcessor(){};

        template<class TTraits, class TDistributionType>
        inline void communicateProcessor(){};

        template<class TTraits>
        inline void runProcessor(int k){};

        inline void setInterfaceID(int id) {mInterfaceID[0]=id;};

        inline void setInterfaceID(const std::vector<int>& id) {mInterfaceID=id;};

        template<class TLattice>
        inline bool apply(int k) {
            if (Geometry<TLattice>::getBoundaryType(k) == -1) return false;
            for (int i : mInterfaceID){
                if(Geometry<TLattice>::getBoundaryType(k) == i) return true;
            }
            return false;
        }

        inline void initialise(Model* model) {mModel = model;};

    private:

        std::vector<int> mInterfaceID = {1};

    protected:
        Model *mModel;
};
