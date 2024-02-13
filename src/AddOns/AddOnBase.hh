#pragma once
#include "../Geometry.hh"

class AddOnBase{
    
    public:

        template<class TTraits>
        inline void compute(int k);

        template<class TTraits>
        inline void communicate();

        inline void setNodeID(int id) {mNodeID = {id};};
        inline void setNodeID(std::vector<int> id) {mNodeID = id;};

        template<class TLattice>
        inline bool apply(int k) {
            if (Geometry<TLattice>::getBoundaryType(k) == -1) return false;
            for (int i : mNodeID){
                if(Geometry<TLattice>::getBoundaryType(k) == i) return true;
            }
            return false;
        }

    private:

        std::vector<int> mNodeID = {};

};

template<class TTraits>
inline void AddOnBase::compute(int k){

}

template<class TTraits>
inline void AddOnBase::communicate() {
    
}
