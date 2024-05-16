#pragma once
#include "../Geometry.hh"
#include "../Service.hh"  //TMP: Default NodeID warning

class Model;

class BoundaryBase {
   public:
    template <class TTraits, class TDistributionType>
    inline void compute(TDistributionType& mDistribution, int k);

    template <class TTraits>
    inline void communicate(){};

    template <class TTraits>
    inline void communicateProcessor(){};

    template <class TTraits, class TDistributionType>
    inline void communicateProcessor(){};

    template <class TTraits>
    inline void runProcessor(int k){};

    // TMP: Default NodeID warning
    inline void setNodeID(int id, bool preset = false) {
        mNodeID = {id};
        preset_warning = preset;
    };
    inline void setNodeID(const std::vector<int>& id, bool preset = false) {
        mNodeID = id;
        preset_warning = preset;
    };

    template <class TLattice>
    inline bool apply(int k) {
        if (Geometry<TLattice>::getBoundaryType(k) == -1) return false;
        for (int i : mNodeID) {
            // TMP: Default NodeID warning
            if (Geometry<TLattice>::getBoundaryType(k) == i) {
                if (preset_warning) {
#pragma omp critical
                    print("\033[31;1mDEPRECATION WARNING\033[0m: Using default NodeID (", i,
                          ") for a boundary. Please explicitly set NodeIDs for all boundaries in use.");
                    preset_warning = false;
                }
                return true;
            }
        }
        return false;
    }

    inline void initialise(Model* model) { mModel = model; };

    std::vector<int> mInterfaceID;

   private:
    std::vector<int> mNodeID = {};
    bool preset_warning = false;  // TMP: Default NodeID warning

   protected:
    Model* mModel;
};
