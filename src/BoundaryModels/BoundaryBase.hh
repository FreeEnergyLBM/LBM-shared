#pragma once
#include "../Geometry.hh"

class BoundaryBase{
    public:

        template<class TTraits>
        inline void precompute(int k){};

        template<class TTraits>
        inline void communicatePrecompute(){};

        template<class TTraits>
        inline void communicate(){};

        template<class TTraits>
        inline void communicatePostProcess(){};

        template<class TTraits, class TDistributionType>
        inline void communicatePrecompute(){};

        template<class TTraits, class TDistributionType>
        inline void communicatePostProcess(){};

        template<class TTraits>
        inline void postprocess(int k){};

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

    private:

        std::vector<int> mInterfaceID = {1};


};