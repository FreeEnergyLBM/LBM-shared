#pragma once
#include "BoundaryBase.hh"
#include "../Geometry.hh"
#include "../LBModels/ModelBase.hh"
#include <map>
#include <algorithm>


class ZouHe : public BoundaryBase {
    public:
        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        inline void setInterfaceID(int id) {mInterfaceID={id};};
        inline void setInterfaceID(const std::vector<int>& id) {mInterfaceID=id;};

        enum BoundaryTypes { DENSITY=0, VELOCITY=1 };
        int boundaryType;

    private:
        template<class TLattice>
        inline void initialiseBoundaryValue(int k);
        std::map<int,double> mBoundaryValues;

        std::vector<int> mInterfaceID = {1};
};


class ZouHeDensity : public ZouHe {
    public:
        ZouHeDensity() { boundaryType = DENSITY; };
};


class ZouHeVelocity : public ZouHe {
    public:
        ZouHeVelocity() { boundaryType = VELOCITY; };
};


inline bool isIn(int value, std::vector<int> array) {
    for (int iTest: array) {
        if (value == iTest) return true;
    }
    return false;
}


template<class TLattice>
inline void getNormal(int k, int& normalDim, int& normalDir) {
    auto normalVec = BoundaryLabels<TLattice::NDIM>::template get<TLattice>(k).NormalDirection;
    normalDim = std::distance(begin(normalVec), std::find_if(begin(normalVec), end(normalVec), [](int ni) {return ni!=0;}));
    normalDir = normalVec[normalDim];
}


template<class TLattice>
inline void ZouHe::initialiseBoundaryValue(int k) {
    if (boundaryType == DENSITY) {
        mBoundaryValues.insert({k, Density<>::template get<TLattice>(k)});
    } else {
        int normalDim, normalDir;
        getNormal<TLattice>(k, normalDim, normalDir);
        double velOut = normalDir * Velocity<>::template get<TLattice,TLattice::NDIM>(k, normalDim);
        mBoundaryValues.insert({k, velOut});
    }
}


template<class TTraits, class TDistributionType>
inline void ZouHe::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;
    if (!isIn(Geometry<Lattice>::getBoundaryType(k), mInterfaceID)) return;
    if (mBoundaryValues.count(k) == 0) initialiseBoundaryValue<Lattice>(k); // Initialise the fixed value if it is not already

    // Get necessary values
    int normalDim, normalDir;
    getNormal<Lattice>(k, normalDim, normalDir);
    auto distrK = distribution.getDistributionPointer(k);
    auto model = static_cast<ModelBase<Lattice,TTraits>*>(mModel);

    // Sum the distribution function groups and get unknowns
    double distOut = 0;
    double distNeutral = 0;
    std::vector<int> unknowns;
    unknowns.reserve(TTraits::Stencil::Q/2);
    for (int iQ=0; iQ<TTraits::Stencil::Q; iQ++) {
        int cNormal = normalDir * TTraits::Stencil::Ci_xyz(normalDim)[iQ];
        if (cNormal == 0) {
            distNeutral += distrK[iQ];
        } else if (cNormal < 0) {
            distOut += distrK[iQ];
        } else {
            unknowns.push_back(iQ);
        }
    }

    // Compute the unknown velocity or density
    double *density = Density<>::template getAddress<Lattice>(k);
    double *velocity = Velocity<>::template getAddress<Lattice,Lattice::NDIM>(k, normalDim);
    if (boundaryType == DENSITY) {
        *density = mBoundaryValues[k];
        *velocity = normalDir * (1 - (distNeutral + 2*distOut) / *density);
    } else {
        double velOut = mBoundaryValues[k];
        *velocity = normalDir * velOut;
        *density = (distNeutral + 2*distOut) / (1 - velOut);
    }

    // Compute the unknown distribution functions
    for (int iQ: unknowns) {
        int iQOpp = distribution.getOpposite(iQ);
        double distrEq = model->computeEquilibrium(k, iQ);
        double distrEqOpp = model->computeEquilibrium(k, iQOpp);
        distrK[iQ] = distrK[iQOpp] + (distrEq - distrEqOpp);
    }

    // TODO: Add parallel velocities, sort out how equilibrium gets calculated
    // Is special treatment needed at the corners?
}
