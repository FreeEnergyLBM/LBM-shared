#pragma once
#include "../Geometry.hh"

class Model;

class BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);

    template <class TTraits>
    inline void communicate(){};

    template <class TTraits, class TDistributionType>
    inline void communicate(TDistributionType& mDistribution){};

    template <class TTraits>
    inline void communicateProcessor(){};

    template <class TTraits, class TDistributionType>
    inline void communicateProcessor(){};

    template <class TTraits>
    inline void runProcessor(int k){};

    inline void setNodeID(int id) { mNodeID = {id}; };
    inline void setNodeID(const std::vector<int>& id) { mNodeID = id; };

    template <class TLattice>
    inline bool apply(int k) {
        if (Geometry<TLattice>::getBoundaryType(k) == -1) return false;
        for (int i : mNodeID) {
            if (Geometry<TLattice>::getBoundaryType(k) == i) {
                return true;
            }
        }
        return false;
    }

    inline void initialise(Model* model) { mModel = model; };

    std::vector<int> mInterfaceID;

   private:
    std::vector<int> mNodeID = {};

   protected:
    Model* mModel;
};
