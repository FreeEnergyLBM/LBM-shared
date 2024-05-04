#pragma once
#include <vector>

/**
 * \brief The base class containing the distribution function data.
 * This will be inherited from in the Data classes (see Data.hh).
 * The Parameter class contains a static vector and functions to return
 */
template <class TStencil>  // Distribution must know information about the stencil as this determines the size
                           // of the vectors and the data layout
struct Distribution_Base {
    using Stencil = TStencil;

    int mQ = TStencil::Q;
    int ma_Opposites[TStencil::Q];  //!< Array containing the opposite indices at each index (rotated by 180 degrees).
    std::vector<double> mv_Distribution;      //!< Vector that will store the distribution information
    std::vector<double> mv_OldDistribution;   //!< Possibly unused vector storing the old distributions (at
                                              //!< t_old=t_new-1) //This is required in the collision step
    std::vector<double> mv_CommDistribution;  //!< Vector that will store the distribution information for communication
    std::vector<double> mv_EquilibriumDistribution;
    std::vector<int>& mv_DistNeighbors;  //!< Reference to vector containing neighbor information

    Distribution_Base(std::vector<int>& neighbors) : mv_DistNeighbors(neighbors) {
        // Calculate the k offset for the neighbors in each direction
        for (int idx = 0; idx < TStencil::Q; idx++) {
            ma_Opposites[idx] = TStencil::Opposites[idx];
        }
    }

    inline void saveEquilibriums(const double* equilibrium, int k) {}

    //! Returns the opposite index at the chosen index (Rotation by 180 degrees).
    inline int getOpposite(int idx) { return ma_Opposites[idx]; }

    //! Get a vector containing the distributions
    inline std::vector<double>& getDistribution() { return mv_Distribution; }

    //! Get a vector containing the distributions for communication (along direction `neighbor`?)
    inline std::vector<double>& getCommDistribution() { return mv_CommDistribution; }

    //! Get a constant pointer to the the distribution at lattice point k and pointing in direction 0
    inline const double* getDistributionPointer(const int k) const { return &mv_Distribution[k * TStencil::Q]; }

    //! Get a pointer to the the distribution at lattice point k and pointing in direction 0
    inline double* getDistributionPointer(const int k) { return &mv_Distribution[k * TStencil::Q]; }

    //! Get const distribution value at a given index
    inline const double& getDistribution(const int idx) const { return mv_Distribution[idx]; }

    //! Get distribution value at a given index
    inline double& getDistribution(const int idx) { return mv_Distribution[idx]; }

    //! Get old distribution value at a given index
    inline std::vector<double>& getDistributionOld() { return mv_OldDistribution; }

    inline const double* getDistributionOldPointer(const int k) const { return &mv_OldDistribution[k * TStencil::Q]; }

    inline double* getDistributionOldPointer(const int k) { return &mv_OldDistribution[k * TStencil::Q]; }

    inline const double& getDistributionOld(const int k) const { return mv_OldDistribution[k]; }

    inline double& getDistributionOld(const int k) { return mv_OldDistribution[k]; }

    //! Get a vector containing the distributions
    inline std::vector<double>& getEquilibrium() { return mv_EquilibriumDistribution; }

    //! Get a constant pointer to the the distribution at lattice point k and pointing in direction 0
    inline const double* getEquilibriumPointer(const int k) const {
        return &mv_EquilibriumDistribution[k * TStencil::Q];
    }

    //! Get a pointer to the the distribution at lattice point k and pointing in direction 0
    inline double* getEquilibriumPointer(const int k) { return &mv_EquilibriumDistribution[k * TStencil::Q]; }

    //! Get const distribution value at a given index
    inline const double& getEquilibrium(const int idx) const { return mv_EquilibriumDistribution[idx]; }

    //! Get distribution value at a given index
    inline double& getEquilibrium(const int idx) { return mv_EquilibriumDistribution[idx]; }
};
