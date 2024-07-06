#pragma once
#include "../Geometry.hh"

class Model;

class AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

    template <class TTraits>
    inline void communicate();

    inline void initialise(Model* model) { mModel = model; };

    inline void setNodeID(int id) { mNodeID = {id}; };

    inline void setNodeID(std::vector<int> id) { mNodeID = id; };

    template <class TLattice>
    inline bool apply(int k) {
        const int& boundarytype = Geometry<TLattice>::getBoundaryType(k);

        if (boundarytype == -1) return false;
        for (int i : mNodeID) {
            if (boundarytype == i) return true;
        }
        return false;
    }

   private:
    std::vector<int> mNodeID = {};

   protected:
    Model* mModel;
};

template <class TTraits>
inline void AddOnBase::compute(int k) {}

template <class TTraits>
inline void AddOnBase::communicate() {}
