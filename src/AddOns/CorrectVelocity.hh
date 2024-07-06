#pragma once

class CorrectVelocity : public AddOnBase {
   public:
    template <class TTraits>
    inline void compute(int k);

   private:
    std::vector<double> mWallVelocity = {0, 0, 0};
    int mSolidPhase = 1;
};

template <class TTraits>
inline void CorrectVelocity::compute(int k) {
    using Lattice = typename TTraits::Lattice;
    const int NDIM = Lattice::NDIM;

    const double& orderparam = getInstance<OrderParameter, 10, Lattice>(mSolidPhase)[k];

    Velocity<>::get<Lattice, NDIM>(k, 0) +=
        1.0 / 2.0 * orderparam * (mWallVelocity[0] - Velocity<>::get<Lattice, NDIM>(k, 0));
    if constexpr (NDIM >= 2) {
        Velocity<>::get<Lattice, NDIM>(k, 1) +=
            1.0 / 2.0 * orderparam * (mWallVelocity[0] - Velocity<>::get<Lattice, NDIM>(k, 1));
    }
    if constexpr (NDIM >= 3) {
        Velocity<>::get<Lattice, NDIM>(k, 2) +=
            1.0 / 2.0 * orderparam * (mWallVelocity[0] - Velocity<>::get<Lattice, NDIM>(k, 2));
    }
}
