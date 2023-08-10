#pragma once

class BoundaryBase{
    public:

        template<class TTraits>
        inline void precompute(int k){};

        template<class TTraits>
        inline void communicatePrecompute(){};

        template<class TTraits>
        inline void communicatePostProcess(){};

        template<class TTraits, class TDistributionType>
        inline void communicatePrecompute(TDistributionType& mDistribution){};

        template<class TTraits, class TDistributionType>
        inline void communicatePostProcess(TDistributionType& mDistribution){};

        template<class TTraits>
        inline void postprocess(int k){};

    private:


};