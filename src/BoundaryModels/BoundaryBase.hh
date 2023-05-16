#pragma once

class BoundaryBase{
    public:

        //virtual void compute(auto& m_Distribution, int k, int idx) const;

        inline virtual void precompute(const int k);

        inline virtual double computeDensitySource(const int k) const;

        inline virtual double computeVelocitySource(const int xyz, const int k) const;

        inline virtual void postprocess(const int k);

    private:


};

inline void BoundaryBase::precompute(const int k) {
    
}

inline void BoundaryBase::postprocess(const int k) {
    
}

inline double BoundaryBase::computeDensitySource(const int k) const {

    return 0;

}

inline double BoundaryBase::computeVelocitySource(const int xyz,const int k) const {

    return 0;

}