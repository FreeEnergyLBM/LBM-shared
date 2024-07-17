#pragma once
#include <dirent.h>

#include <functional>
#include <map>

#include "Data.hh"
#include "Parameters.hh"
#include "Stencil.hh"

/**
 * \file Geometry.hh
 * \brief This is not fully implemented but will contain information relevant to boundary conditions
 * This contains a geometry class that will contain information about the solid geometry and boundary locations
 * in the simulation.
 */

/**
 * \brief Geometry contains functions to initialise and access the chosen geometry.
 * This class contains functions to initialise the chosen geometry and determine whether the current TLattice
 * point is a solid or a periodic boundary.
 * \tparam Empty template parameter, used to ensure GETPROPERTIES is defined before it is needed in this function.
 */
template <class TLattice>
class Geometry {
   public:
    /**
     * \brief Returns true if the current TLattice point lies on a solid boundary.
     * \param k Index of current TLattice point.
     * \return True if TLattice point lies on a solid
     */
    static inline bool isBoundary(int k);

    /**
     * \brief Returns the normal to the boundary as a cartesian vector.
     * \tparam TStencil The type of the velocity stencil.
     * \param condition Function that returns the boundary index for a given lattice index k.
     * \param neighbors Vector that contains the neighboring index in each velocity direction on each lattice index.
     * \param fluidvals Boundary index of lattice points that are fluids (for normal calculation).
     * \param k Latice index.
     */
    template <class TStencil>
    static inline std::array<int8_t, TLattice::NDIM> findNormal(std::function<int(const int)> condition,
                                                                const std::vector<int>& neighbors,
                                                                const std::vector<int> fluidvals, int k);

    /**
     * \brief Returns true if the current lattice point is a corner boundary.
     * \param condition Function that returns the boundary index for a given lattice index k.
     * \param normal Normal to the boundary as a cartesian vector.
     * \param fluidvals Boundary index of lattice points that are fluids (for normal calculation).
     * \param k Latice index.
     */
    static inline bool isCorner(std::function<int(const int)> condition,
                                const std::array<int8_t, TLattice::NDIM>& normal, const std::vector<int> fluidvals,
                                int k);

    /**
     * \brief Returns true if the current lattice point is a corner boundary.
     * \param k Latice index.
     */
    static inline bool isCorner(int k);

    /**
     * \brief Returns true if the current lattice point is surrounded by boundary/other bulk solid nodes..
     * \param condition Function that returns the boundary index for a given lattice index k.
     * \param normal Normal to the boundary as a cartesian vector.
     * \param fluidvals Boundary index of lattice points that are fluids (for normal calculation).
     * \param k Latice index.
     */
    static inline bool isBulkSolid(std::function<int(const int)> condition, const std::vector<int>& neighbors,
                                   std::vector<int> fluidvals, int k);

    /**
     * \brief Returns true if the current lattice point is surrounded by boundary/other bulk solid nides.
     * \param k Latice index.
     */
    static inline bool isBulkSolid(int k);

    /**
     * \brief Initialises the BoundaryLabels class and sets fluid boundary indices.
     * \param condition Function that returns the boundary index for a given lattice index k.
     * \param fluidval Boundary index of lattice points that are fluids (for normal calculation).
     */
    static inline void initialiseBoundaries(std::function<int(const int)> condition, int fluidval = 0);

    /**
     * \brief Initialises the BoundaryLabels class and sets fluid boundary indices.
     * \param condition Function that returns the boundary index for a given lattice index k.
     * \param fluidval Boundary indices of lattice points that are fluids (for normal calculation).
     */
    static inline void initialiseBoundaries(std::function<int(const int)> condition, std::vector<int> fluidvals);

    /**
     * \brief Initialises the BoundaryLabels class and sets fluid boundary indices.
     * \param indices Vector that contains the boundary indices for the entire domain.
     * \param fluidvals Boundary indices of lattice points that are fluids (for normal calculation).
     */
    static inline void initialiseBoundaries(std::vector<int> indices, std::vector<int> fluidvals);

    /**
     * \brief Initialises the BoundaryLabels class and sets fluid boundary indices.
     * \param indices Vector that contains the boundary indices for the entire domain.
     * \param fluidval (optional) Boundary index of lattice points that are fluids (for normal calculation).
     */
    static inline void initialiseBoundaries(std::vector<int> indices, int fluidval = 0);

    /**
     * \brief Returns the boundary index at the given lattice index.
     * \param k Lattice index.
     */
    static inline int& getBoundaryType(int k);

    // New function declaration to modify the boundary label
    static inline void modifyBoundaryLabel(int k, int modifiedLabel);

    static std::vector<int> mFluidVals;  //!< Boundary indices corresponding to fluid nodes.

   private:
    enum {
        BulkSolid = -1,
        Fluid = 0,
        Wall = 1,
        Wall2 = 2,
        InletWall = 3,
        OutletWall = 4,
        HumidityInterface = 5,
        RefillNode = 6
    };

    /**
     * \brief Returns true if the lattice index is not a fluid node.
     * \param condition Function that returns the boundary index for a given lattice index k.
     * \param fluidvals Boundary index of lattice points that are fluids (for normal calculation).
     * \param k Latice index.
     */
    static inline bool apply(std::function<int(const int)> condition, std::vector<int> fluidvals, int k) {
        for (int i : fluidvals) {
            if (condition(k) == i) return false;
        }
        return true;
    }
};

template <class TLattice>
std::vector<int> Geometry<TLattice>::mFluidVals;

/**
 * \details This function iterates through all velocity directions and then iterates through the fluid boundary indices.
 *          If it finds that the boundary indiex in any velocity direction is a fluid, it returns false. Otherwise it
 * returns true
 */
template <class TLattice>
inline bool Geometry<TLattice>::isBulkSolid(std::function<int(const int)> condition, const std::vector<int>& neighbors,
                                            std::vector<int> fluidvals, int k) {
    mFluidVals = fluidvals;

    using Stencil = std::conditional_t<TLattice::NDIM == 1, D1Q3, std::conditional_t<TLattice::NDIM == 2, D2Q9, D3Q27>>;

    bool bulksolid = true;

    for (int idx = 0; idx < Stencil::Q; idx++) {
        for (int i : fluidvals) {
            if (  // condition(neighbors[k * Stencil::Q + idx]) != condition(k) ||
                condition(neighbors[k * Stencil::Q + idx]) == i)
                bulksolid = false;
        }
    }

    return bulksolid;
}

/**
 * \details This function iterates through all lattice points in the simulation domain and initialises the
 *          array of structures stored in the BoundaryLabels class with the boundary label, whether or not the
 *          point is a corner and the normal direction pointing into the fluid.
 */
template <class TLattice>
inline void Geometry<TLattice>::initialiseBoundaries(std::function<int(const int)> condition,
                                                     const std::vector<int> fluidvals) {
    using Stencil = std::conditional_t<TLattice::NDIM == 1, D1Q3, std::conditional_t<TLattice::NDIM == 2, D2Q9, D3Q27>>;

    using data = Data_Base<TLattice, Stencil>;
    data::getInstance().generateNeighbors();
    std::vector<int>& neighbors = data::getInstance().getNeighbors();

    for (int k : RangeK<TLattice>()) {
        int solidval;
        if (isBulkSolid(condition, neighbors, fluidvals, k)) {
            solidval = -1;
        } else {
            solidval = condition(k);
        }

        std::array<int8_t, TLattice::NDIM> normal = findNormal<Stencil>(condition, neighbors, fluidvals, k);
        Boundary<TLattice::NDIM> boundaryk = {solidval, isCorner(condition, normal, fluidvals, k), normal};
        BoundaryLabels<TLattice::NDIM>::template initialise<TLattice>(boundaryk, k);
    }

    neighbors.resize(0);
}

/**
 * \details This class calls the InitialiseBoundaries class after converting the integer fluidval to a vector.
 */
template <class TLattice>
inline void Geometry<TLattice>::initialiseBoundaries(std::function<int(const int)> condition, int fluidval) {
    initialiseBoundaries(condition, std::vector<int>{fluidval});
}

/**
 * \details This function first converts the index array to a function before calling the initialiseBoundaries
 * implementation.
 */
template <class TLattice>
inline void Geometry<TLattice>::initialiseBoundaries(std::vector<int> indices, const std::vector<int> fluidvals) {
    auto indexFunction = [&indices](const int k) { return indices[computeKGlobal<TLattice>(k)]; };
    initialiseBoundaries(indexFunction, fluidvals);
}

template <class TLattice>
inline void Geometry<TLattice>::initialiseBoundaries(std::vector<int> indices, int fluidval) {
    initialiseBoundaries(indices, std::vector<int>{fluidval});
}

/**
 * \details This function evaluates an if statement to determine if we are on the solid boundary. This boundary
 *          is currently determined by this function, but in the future it will be chosen elsewhere so that no
 *          modification to the source file is needed. Returns true if current TLattice point is a solid.
 */
template <class TLattice>
inline bool Geometry<TLattice>::isBoundary(int k) {
    for (int i : mFluidVals) {
        if (BoundaryLabels<TLattice::NDIM>::template get<TLattice>(k).Id == i) return false;
    }
    return true;
}

/**
 * \details This function returns the true/false value stored in BoundaryLabels indicating whether the current
 *          lattice point is a corner boundary,
 */
template <class TLattice>
inline bool Geometry<TLattice>::isCorner(int k) {
    return BoundaryLabels<TLattice::NDIM>::template get<TLattice>(k).IsCorner;
}

/**
 * \details This function returns true if the boundary index stored in BoundaryLabels at the given lattice
 *          index is equal to -1.
 */
template <class TLattice>
inline bool Geometry<TLattice>::isBulkSolid(int k) {
    return (BoundaryLabels<TLattice::NDIM>::template get<TLattice>(k).Id == -1);
}

/**
 * \details This function finds the normal by computing the sum of the number of fluid nodes adjacent in each
 *          cartesian direction. The normal in that direction is then -1, 0 or +1 depending on whether the sum is <0,
 * ==0 or <0. It first sums the indices around the the current nodes and incriments the sum in each direction if the
 * neighboring node is a fluid and is not the same boundary index as the current node. It then checks if there is a
 * difference in the sum in each direction, If the normal would not form a 45 or 90 degree angle, it sets the sum to
 * zero in the direction where the sum is smallest. It then sets the normal to the sign of the sum in each direction.
 */
template <class TLattice>
template <class TStencil>
inline std::array<int8_t, TLattice::NDIM> Geometry<TLattice>::findNormal(std::function<int(const int)> condition,
                                                                         const std::vector<int>& neighbors,
                                                                         const std::vector<int> fluidvals, int k) {
    std::array<int8_t, TLattice::NDIM> normal = {};

    std::vector<int> sum(TLattice::NDIM, 0);

    for (int idx = 0; idx < TStencil::Q; idx++) {
        if (!apply(condition, fluidvals, neighbors[k * TStencil::Q + idx]) &&
            condition(k) != condition(neighbors[k * TStencil::Q + idx])) {
            sum[0] += TStencil::Ci_xyz(0)[idx];
            if constexpr (TLattice::NDIM > 1) sum[1] += TStencil::Ci_xyz(1)[idx];
            if constexpr (TLattice::NDIM > 2) sum[2] += TStencil::Ci_xyz(2)[idx];
        }
    }

    if constexpr (TLattice::NDIM > 1) {
        if (abs(sum[0]) > abs(sum[1])) {
            if constexpr (TLattice::NDIM > 2) {
                if (abs(sum[0]) > abs(sum[2])) {
                    sum[2] = 0;
                } else if (abs(sum[0]) < abs(sum[2])) {
                    sum[0] = 0;
                }
                if (abs(sum[1]) > abs(sum[2])) {
                    sum[2] = 0;
                } else if (abs(sum[1]) < abs(sum[2])) {
                    sum[1] = 0;
                }
            }
            sum[1] = 0;
        } else if (abs(sum[0]) < abs(sum[1])) {
            if constexpr (TLattice::NDIM > 2) {
                if (abs(sum[0]) > abs(sum[2])) {
                    sum[2] = 0;
                } else if (abs(sum[0]) < abs(sum[2])) {
                    sum[0] = 0;
                }
                if (abs(sum[1]) > abs(sum[2])) {
                    sum[2] = 0;
                } else if (abs(sum[1]) < abs(sum[2])) {
                    sum[1] = 0;
                }
            }
            sum[0] = 0;
        }
    }

    for (int xyz = 0; xyz < TLattice::NDIM; xyz++) {
        normal[xyz] = (int8_t)((sum[xyz] == 0) ? 0 : ((sum[xyz] > 0) ? 1 : -1));
    }

    return normal;
}

/**
 * \details This function returns false if the current lattice node is a fluid. It returns true if more than
 *          one of the normal directions is not equal to zero. Otherwise it returns false.
 */
template <class TLattice>
inline bool Geometry<TLattice>::isCorner(std::function<int(const int)> condition,
                                         const std::array<int8_t, TLattice::NDIM>& normal,
                                         const std::vector<int> fluidvals, int k) {  // DOesnt

    if (!apply(condition, fluidvals, k)) return false;

    int normalsum = 0;

    for (int i = 0; i < (int)normal.size(); i++) {
        normalsum += abs((int)normal[i]);
    }

    return normalsum > 1;
}

/**
 * \details This function returns the boundary index stored in BoundaryLabels at the given lattice index.
 */
template <class TLattice>
inline int& Geometry<TLattice>::getBoundaryType(int k) {
    return BoundaryLabels<TLattice::NDIM>::template get<TLattice>(k).Id;
}

template <class TLattice>
inline void Geometry<TLattice>::modifyBoundaryLabel(int k, int modifiedLabel) {
    BoundaryLabels<TLattice::NDIM>::template get<TLattice>(k).Id = modifiedLabel;
}
