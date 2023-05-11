#pragma once

template<typename placeholder = void>
class BoundaryBaseTemplate{
    public:

        //virtual void compute(auto& m_Distribution, int k, int idx) const;

        inline virtual void precompute(const int k);

        inline virtual double computeDensitySource(const int k) const;

        inline virtual double computeVelocitySource(const int xyz, const int k) const;

        inline virtual void postprocess(const int k);

    private:


};

template<typename placeholder>
inline void BoundaryBaseTemplate<placeholder>::precompute(const int k) {
    
}

template<typename placeholder>
inline void BoundaryBaseTemplate<placeholder>::postprocess(const int k) {
    
}

template<typename placeholder>
inline double BoundaryBaseTemplate<placeholder>::computeDensitySource(const int k) const {

    return 0;

}

template<typename placeholder>
inline double BoundaryBaseTemplate<placeholder>::computeVelocitySource(const int xyz,const int k) const {

    return 0;

}

typedef BoundaryBaseTemplate<> BoundaryBase;