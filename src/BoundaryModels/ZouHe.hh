#pragma once
#include "BoundaryBase.hh"
#include "../Geometry.hh"
#include "../LBModels/ModelBase.hh"
#include <map>
#include <algorithm>


//ZouHe.hh: Contains the Zou-He boundary condition for setting a constant density or velocity (https://doi.org/10.1063/1.869307).
//Use ZouHeDensity for constant density boundaries and ZouHeVelocity for constant velocity.


class ZouHe : public BoundaryBase {
    public:
        ZouHe() : BoundaryBase(4) {}
        template<class TTraits, class TDistributionType>
        inline void compute(TDistributionType& mDistribution, int k);

        // Used to set a non-perpendicular velocity for constant density boundaries
        void setAngledVelocity(std::vector<double> velocityDirection);

        enum BoundaryTypes { DENSITY=0, VELOCITY=1 };
        int boundaryType;

    private:
        template<class TLattice>
        inline void initialiseBoundaryValue(int k);
        std::map<int,double> mBoundaryValues;

        template<class TLattice>
        inline void angleVelocity(int k, int normalDim);
        bool mUseAngledVelocity = false;
        std::vector<double> mVelocityDirection;
};


class ZouHeDensity : public ZouHe {
    public:
        ZouHeDensity() { boundaryType = DENSITY; };
};


class ZouHeVelocity : public ZouHe {
    public:
        ZouHeVelocity() { boundaryType = VELOCITY; };
};


template<class TLattice>
inline void getNormal(int k, int& normalDim, int& normalDir) {
    auto normalVec = BoundaryLabels<TLattice::NDIM>::template get<TLattice>(k).NormalDirection;
    normalDim = std::distance(begin(normalVec), std::find_if(begin(normalVec), end(normalVec), [](int ni) {return ni!=0;}));
    normalDir = normalVec[normalDim];
}


template<class TLattice>
inline void ZouHe::initialiseBoundaryValue(int k) {
    if (boundaryType == DENSITY) {
        #pragma omp critical
        mBoundaryValues.insert({k, Density<>::template get<TLattice>(k)});
    } else {
        int normalDim, normalDir;
        getNormal<TLattice>(k, normalDim, normalDir);
        double velOut = normalDir * Velocity<>::template get<TLattice,TLattice::NDIM>(k, normalDim);
        #pragma omp critical
        mBoundaryValues.insert({k, velOut});
    }
}


template<class TTraits, class TDistributionType>
inline void ZouHe::compute(TDistributionType& distribution, int k) {
    using Lattice = typename TTraits::Lattice;
    if (!apply<Lattice>(k)) return;
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
        if (mUseAngledVelocity) angleVelocity<Lattice>(k, normalDim);
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

    // TODO: Pressure models
}


void ZouHe::setAngledVelocity(std::vector<double> velocityDirection) {
    mUseAngledVelocity = true;
    mVelocityDirection = velocityDirection;
}


template<class TLattice>
inline void ZouHe::angleVelocity(int k, int normalDim) {
    double *velocity = Velocity<>::template getAddress<TLattice,TLattice::NDIM>(k, 0);
    double vPerp = velocity[normalDim];
    // (De)project the perpendicular velocity onto the velocity direction
    double projectionFactor = 1 / mVelocityDirection[normalDim];
    for (int iDim=0; iDim<TLattice::NDIM; iDim++) {
        velocity[iDim] = vPerp * mVelocityDirection[iDim] * projectionFactor;
    }
}
