#pragma once

class BoundaryBase{
    public:

        template<class T_traits>
        inline void precompute(int k){};

        template<class T_traits>
        inline void communicatePrecompute(){};

        template<class T_traits>
        inline void communicatePostProcess(){};

        template<class T_traits>
        inline void postprocess(int k){};

    private:


};