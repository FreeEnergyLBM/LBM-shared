#pragma once
#include "../Geometry.hh"
#include "../Template.hh"

template <template <class> class TGradientType, class TDirections = Cartesian>
struct GradientBase {
    template <class TTraits, class TParameter>
    inline double compute(const int direction, const int k);

    template <class TObj>
    using GradientType = TGradientType<TObj>;

    template <class TStencil>
    inline static constexpr int getNumberOfDirections() {
        if constexpr (std::is_same_v<TDirections, Cartesian>)
            return TStencil::D;
        else if constexpr (std::is_same_v<TDirections, AllDirections>)
            return TStencil::Q;
        else if constexpr (std::is_same_v<TDirections, One>)
            return 1;
        else
            return TStencil::D;
    }

    template <class TLattice>
    inline bool isBoundary(int k) {
        //if (Geometry<TLattice>::getBoundaryType(k) == -1) return true;
        for (int i : mBoundaryID) {
            // TMP: Default BoundaryID warning
            if (Geometry<TLattice>::getBoundaryType(k) == i) {
                if (preset_warning) {
#pragma omp critical
                    print(
                        "\033[31;1mDEPRECATION WARNING - FIX IMMEDIATELY OR SCRIPT WILL BREAK! \033[0m: "
                        "Using default BoundaryID (",
                        i,
                        ") for a gradient. Please explicitly set BoundaryIDs for all gradients that check boundaries.");
                    preset_warning = false;
                }
                return true;
            }
        }
        return false;
    }

    inline void setBoundaryID(int id, bool preset = false) {
        mBoundaryID = {id};
        preset_warning = preset;
    };
    inline void setBoundaryID(const std::vector<int>& id, bool preset = false) {
        mBoundaryID = id;
        preset_warning = preset;
    };

    std::vector<int> mBoundaryID = {1};
    bool preset_warning = true;

    inline void setPrefactor(double prefactor) { mPrefactor = prefactor; }

    double mPrefactor = 0;

    inline void setInterfaceDistance(double (*distance)(int k, int idx)) {}

    inline void setInterfaceVal(double value) {}
};
