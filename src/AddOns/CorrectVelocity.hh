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
    Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0) +=
        1.0 / 2.0 *
        OrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice>(k, mSolidPhase) *
        (mWallVelocity[0] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 0));
    if constexpr (TTraits::Lattice::NDIM >= 2)
        Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 1) +=
            1.0 / 2.0 *
            OrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice>(k, mSolidPhase) *
            (mWallVelocity[0] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 1));
    if constexpr (TTraits::Lattice::NDIM >= 3)
        Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 2) +=
            1.0 / 2.0 *
            OrderParameter<TTraits::NumberOfComponents - 1>::template get<typename TTraits::Lattice>(k, mSolidPhase) *
            (mWallVelocity[0] - Velocity<>::get<typename TTraits::Lattice, TTraits::Lattice::NDIM>(k, 2));
}