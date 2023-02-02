#include <vector>
using namespace std;
template<typename T>
struct D2Q9{
    static const int m_D=2;
    static const int m_Q=9;
    static constexpr double m_Cs2=0.57735026919;
    static constexpr int ma_Cx[m_Q]={0,1,-1,0,0,1,-1,1,-1};
    static constexpr int ma_Cy[m_Q]={0,0,0,1,-1,1,-1,-1,1};
    static constexpr int ma_Cz[m_Q]={0,0,0,0,0,0,0,0,0};
    template<int d>
    static constexpr auto ma_Cxyz()->const int(&)[m_Q]{
        if constexpr (d==0) {
            return ma_Cx;
        }
        if constexpr (d==1) {
            return ma_Cy;
        }
        if constexpr (d==2) {
            return ma_Cz;
        }
    }
    static constexpr int ma_Opposites[m_Q]={0,2,1,4,3,6,5,8,7};
    static constexpr T ma_Weights[m_Q]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
    template<int idx>
    static const int m_CModulus=ma_Cx[idx]*ma_Cx[idx]+ma_Cy[idx]*ma_Cy[idx];
    //static constexpr int ma_CMod[m_Q]={m_CMod<0>,m_CMod<1>,m_CMod<2>,m_CMod<3>,m_CMod<4>,m_CMod<5>,m_CMod<6>,m_CMod<7>,m_CMod<8>,m_CMod<9>};

    template<int idx>
    static constexpr int ma_Moments[m_Q]={1,
                                    3*(m_CModulus<idx>)-4,
                                    (9*(m_CModulus<idx>*m_CModulus<idx>)-21*(m_CModulus<idx>)+8)/2,
                                    ma_Cx[idx],
                                    (3*(m_CModulus<idx>)-5)*ma_Cx[idx],
                                    ma_Cy[idx],
                                    (3*(m_CModulus<idx>)-5)*ma_Cy[idx],
                                    ma_Cx[idx]*ma_Cx[idx]-ma_Cy[idx]*ma_Cy[idx],
                                    ma_Cx[idx]*ma_Cy[idx]};
    inline static vector<int> mv_MRTWeights(const T& invtau){
        return {0,1,1,0,1,0,1,invtau,invtau};
    }
};

template<typename T>
struct D2Q5{
    static const int m_D=2;
    static const int m_Q=5;
    static constexpr double m_Cs2=0.57735026919;
    static constexpr int ma_Cx[m_Q]={0,1,-1,0,0};
    static constexpr int ma_Cy[m_Q]={0,0,0,1,-1};
    static constexpr int ma_Cz[m_Q]={0,0,0,0,0};
    static constexpr int ma_Opposites[m_Q]={0,2,1,4,3};
    static constexpr T ma_Weights[m_Q]={1/3.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0};
    template<int idx>
    static const int m_CModulus=ma_Cx[idx]*ma_Cx[idx]+ma_Cy[idx]*ma_Cy[idx];
    //static constexpr int ma_CMod[m_Q]={m_CMod<0>,m_CMod<1>,m_CMod<2>,m_CMod<3>,m_CMod<4>,m_CMod<5>,m_CMod<6>,m_CMod<7>,m_CMod<8>,m_CMod<9>};

    template<int idx>
    static constexpr int ma_Moments[m_Q]={1,
                                    3*(m_CModulus<idx>)-4,
                                    ma_Cx[idx],
                                    ma_Cy[idx],
                                    ma_Cx[idx]*ma_Cx[idx]-ma_Cy[idx]*ma_Cy[idx]};
    inline static vector<int> mv_MRTWeights(const T& invtau){
        return {0,0,0,1,invtau};
    }
};

template<typename T>
struct D3Q19{
    static constexpr int m_D=3;
    static constexpr int m_Q=19;
    static constexpr double m_Cs2=0.57735026919;
    static constexpr int ma_Cx[m_Q]={0,1,-1,0,0,0,0,1,-1,1,-1,0,0,0,0,1,-1,1,-1};
    static constexpr int ma_Cy[m_Q]={0,0,0,1,-1,0,0,1,-1,-1,1,1,-1,1,-1,0,0,0,0};
    static constexpr int ma_Cz[m_Q]={0,0,0,0,0,1,-1,0,0,0,0,1,-1,-1,1,1,-1,-1,1};
    static constexpr int ma_Opposites[m_Q]={0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17};
    static constexpr T ma_Weights[m_Q]={1.0/3.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
				                        1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/36.0,
				                        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
				                        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
				                        1.0/36.0, 1.0/36.0, 1.0/36.0};
    template<int idx>
    static const int m_CModulus=ma_Cx[idx]*ma_Cx[idx]+ma_Cy[idx]*ma_Cy[idx]+ma_Cz[idx]*ma_Cz[idx];
    //static constexpr int ma_CMod[m_Q]={m_CMod<0>,m_CMod<1>,m_CMod<2>,m_CMod<3>,m_CMod<4>,m_CMod<5>,m_CMod<6>,m_CMod<7>,m_CMod<8>,m_CMod<9>};

    template<int idx>
    static constexpr int ma_Moments[m_Q]={1,
                                    19*(m_CModulus<idx>)-30,
                                    (21*(m_CModulus<idx>*m_CModulus<idx>)-53*(m_CModulus<idx>)+24)/2,
                                    ma_Cx[idx],
                                    (5*(m_CModulus<idx>)-9)*ma_Cx[idx],
                                    ma_Cy[idx],
                                    (5*(m_CModulus<idx>)-9)*ma_Cy[idx],
                                    ma_Cz[idx],
                                    (5*(m_CModulus<idx>)-9)*ma_Cz[idx],
                                    3*ma_Cx[idx]*ma_Cx[idx]-m_CModulus<idx>,
                                    (3*m_CModulus<idx>-5)*(3*ma_Cx[idx]*ma_Cx[idx]-m_CModulus<idx>),
                                    ma_Cy[idx]*ma_Cy[idx]-ma_Cz[idx]*ma_Cz[idx],
                                    (3*m_CModulus<idx>-5)*(ma_Cy[idx]*ma_Cy[idx]-ma_Cz[idx]*ma_Cz[idx]),
                                    ma_Cx[idx]*ma_Cy[idx],
                                    ma_Cy[idx]*ma_Cz[idx],
                                    ma_Cx[idx]*ma_Cz[idx],
                                    ma_Cx[idx]*(ma_Cy[idx]*ma_Cy[idx]-ma_Cz[idx]*ma_Cz[idx]),
                                    ma_Cy[idx]*(ma_Cz[idx]*ma_Cz[idx]-ma_Cx[idx]*ma_Cx[idx]),
                                    ma_Cz[idx]*(ma_Cx[idx]*ma_Cx[idx]-ma_Cy[idx]*ma_Cy[idx])};
    inline static vector<int> mv_MRTWeights(const T& invtau){
        return {0,1,1,0,1,0,1,0,1,invtau,1,invtau,1,invtau,invtau,invtau,1,1,1};
    }
};

template<typename T>
struct D3Q7{
    static const int m_D=3;
    static const int m_Q=7;
    static constexpr double m_Cs2=1.0/3.0;
    static constexpr int ma_Cx[m_Q]={0,1,-1,0,0,0,0};
    static constexpr int ma_Cy[m_Q]={0,0,0,1,-1,0,0};
    static constexpr int ma_Cz[m_Q]={0,0,0,0,0,1,-1};
    static constexpr int ma_Opposites[m_Q]={0,2,1,4,3,6,5};
    static constexpr T ma_Weights[m_Q]={1.0/4.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0};
    template<int idx>
    static const int m_CModulus=ma_Cx[idx]*ma_Cx[idx]+ma_Cy[idx]*ma_Cy[idx]+ma_Cz[idx]*ma_Cz[idx];
    //static constexpr int ma_CMod[m_Q]={m_CMod<0>,m_CMod<1>,m_CMod<2>,m_CMod<3>,m_CMod<4>,m_CMod<5>,m_CMod<6>,m_CMod<7>,m_CMod<8>,m_CMod<9>};

    template<int idx>
    static constexpr int ma_Moments[m_Q]={1,
                                    ma_Cx[idx],
                                    ma_Cy[idx],
                                    ma_Cz[idx],
                                    6-7*m_CModulus<idx>,
                                    3*ma_Cx[idx]-(m_CModulus<idx>*m_CModulus<idx>),
                                    ma_Cy[idx]*ma_Cy[idx]-ma_Cz[idx]*ma_Cz[idx]};
    inline static vector<int> mv_MRTWeights(const T& invtau){
        return {0,0,0,0,1,1,invtau};
    }
};