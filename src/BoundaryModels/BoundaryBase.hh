#pragma once

class BoundaryBase{
    public:

        template<class traits>
        inline void precompute(const int k){};

        template<class traits>
        inline void communicatePrecompute(){};

        template<class traits>
        inline void communicatePostProcess(){};

        template<class traits>
        inline void postprocess(const int k){};

    private:


};