#pragma once

#include <array>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Global.hh"
#include "Mpi.hh"
#include "Stencil.hh"

// Service.hh: This will contain some commonly used functions with various uses.

// is nan

// stat printing

template <class TLattice>
inline int computeXGlobal(const int k)  // Compute X direction from a given k, the convention in this code is that
                                        // k will iterate over the z direction first, then increment y by 1 once it
                                        // reaches LZ, then repeat the iteration over z. Once it reaches LY x will be
                                        // incremented and this process continues
{
    int x = TLattice::LXMPIOffset + int((k - TLattice::HaloSize) / (float)(TLattice::LZ * TLattice::LY));
    if (k < TLattice::HaloSize) {
        if (x - 1 < 0) return TLattice::LX + x - 1;
        return x - 1;
    }
    return x;
}

inline int computeX(const int& LY, const int& LZ,
                    const int k)  // Compute X direction from a given k, the convention in this code is that
                                  // k will iterate over the z direction first, then increment y by 1 once it reaches
                                  // LZ, then repeat the iteration over z. Once it reaches LY x will be incremented and
                                  // this process continues
{
    return int(k / (float)(LZ * LY));
}

inline int computeY(const int& LY, const int& LZ,
                    const int k)  // Compute Y direction from a given k, this uses information from the X direction
{
    return int((k - computeX(LY, LZ, k) * LZ * LY) / (float)LZ);
}

inline int computeZ(
    const int& LY, const int& LZ,
    const int k)  // Compute Y direction from a given k, this uses information from the X and Y directions
{
    return k - computeX(LY, LZ, k) * LZ * LY - computeY(LY, LZ, k) * LZ;
}

/// Compute the global x coordinate using the local index
template <class TLattice>
inline int computeX(int k) {
    int xlocal = k / (TLattice::LYdiv * TLattice::LZdiv);
    return TLattice::LXMPIOffset + xlocal;
}

/// Compute the global y coordinate using the local index
template <class TLattice>
inline int computeY(int k) {
    int lyz = TLattice::LYdiv * TLattice::LZdiv;
    int xlocal = k / lyz;
    int ylocal = (k - xlocal * lyz) / TLattice::LZdiv;
    return TLattice::LYMPIOffset + ylocal;
}

/// Compute the global z coordinate using the local index
template <class TLattice>
inline int computeZ(int k) {
    int zlocal = k % TLattice::LZdiv;
    return TLattice::LZMPIOffset + zlocal;
}

/// Compute the global x,y,z coordinates using the local index
template <class TLattice>
inline std::array<int, 3> computeXYZ(int k) {
    int lyz = TLattice::LYdiv * TLattice::LZdiv;
    int xlocal = k / lyz;
    int ylocal = (k - xlocal * lyz) / TLattice::LZdiv;
    int zlocal = k % TLattice::LZdiv;
    int xglobal = (xlocal + TLattice::LXMPIOffset - TLattice::HaloXWidth + TLattice::LX) % TLattice::LX;
    int yglobal = (ylocal + TLattice::LYMPIOffset - TLattice::HaloYWidth + TLattice::LY) % TLattice::LY;
    int zglobal = (zlocal + TLattice::LZMPIOffset - TLattice::HaloZWidth + TLattice::LZ) % TLattice::LZ;
    return {xglobal, yglobal, zglobal};
}

/// Compute the local index from the local x,y,z coordinates
template <class TLattice>
inline int computeK(int xLocal, int yLocal, int zLocal) {
    int xProc = xLocal + TLattice::HaloXWidth;
    int yProc = yLocal + TLattice::HaloYWidth;
    int zProc = zLocal + TLattice::HaloZWidth;
    int lyProc = TLattice::LYdiv;
    int lzProc = TLattice::LZdiv;
    return xProc * lyProc * lzProc + yProc * lzProc + zProc;
}

/// Compute the local index from the global x,y,z coordinates (-1 if not local)
template <class TLattice>
inline int computeKFromGlobal(int x, int y, int z) {
    int lxProc = TLattice::LXdiv;
    int lyProc = TLattice::LYdiv;
    int lzProc = TLattice::LZdiv;
    int xProc = x - TLattice::LXMPIOffset + TLattice::HaloXWidth;
    int yProc = y - TLattice::LYMPIOffset + TLattice::HaloYWidth;
    int zProc = z - TLattice::LZMPIOffset + TLattice::HaloZWidth;
    // Ensure periodic boundaries are correct and check if not local
    if (xProc < 0) xProc += TLattice::LX;
    if (yProc < 0) yProc += TLattice::LY;
    if (zProc < 0) zProc += TLattice::LZ;
    if (xProc >= lxProc || yProc >= lyProc || zProc >= lzProc) return -1;
    return xProc * lyProc * lzProc + yProc * lzProc + zProc;
}

/// Compute the global index from the local index
template <class TLattice>
inline int computeKGlobal(int k) {
    std::array<int, 3> xyz = computeXYZ<TLattice>(k);
    return xyz[0] * TLattice::LY * TLattice::LZ + xyz[1] * TLattice::LZ + xyz[2];
}

template <class TLattice, typename TOut>
struct RangeXYZIterator {
    RangeXYZIterator(int k) : k(k) {
        auto xyz = computeXYZ<TLattice>(k);
        x = xyz[0];
        y = xyz[1];
        z = xyz[2];
    }

    TOut operator*() const {
        if constexpr (std::is_same<TOut, int>::value) {
            return k;
        } else if constexpr (std::is_same<TOut, std::array<int, 3>>::value) {
            return {x, y, z};
        } else if constexpr (std::is_same<TOut, std::array<int, 4>>::value) {
            return {x, y, z, k};
        }
    }

    // Increment x, y, z, and k, skipping the halo nodes
    RangeXYZIterator<TLattice, TOut> operator++() {
        if (z < zEnd - 1) {
            z++;
            k++;
        } else if (y < yEnd - 1) {
            z = zStart;
            y++;
            k += 1 + 2 * TLattice::HaloZWidth;
        } else if (x < xEnd - 1) {
            z = zStart;
            y = yStart;
            x++;
            k += 1 + 2 * (TLattice::HaloZWidth + TLattice::HaloYWidth * TLattice::LZdiv);
        } else {
            k++;  // Increment to match end() value
        }
        return *this;
    };

    bool operator!=(const RangeXYZIterator& other) const { return k != other.k; };

   private:
    int k, x, y, z;
    int xStart = TLattice::LXMPIOffset;
    int yStart = TLattice::LYMPIOffset;
    int zStart = TLattice::LZMPIOffset;
    int xEnd = xStart + TLattice::subArray[0];
    int yEnd = yStart + TLattice::subArray[1];
    int zEnd = zStart + TLattice::subArray[2];
};

/**
 * \brief Iterator over the index k of the local lattice.
 */
template <class TLattice>
class RangeK {
   public:
    RangeK<TLattice>() {
        kStart = computeK<TLattice>(0, 0, 0);
        kEnd = computeK<TLattice>(TLattice::subArray[0] - 1, TLattice::subArray[1] - 1, TLattice::subArray[2] - 1) + 1;
    }

    using Iterator = RangeXYZIterator<TLattice, int>;
    Iterator begin() { return Iterator(kStart); }
    Iterator end() { return Iterator(kEnd); }

   private:
    int kStart, kEnd;
};

/**
 * \brief Iterator over the x,y,z coordinates of the local lattice.
 */
template <class TLattice>
class RangeXYZ {
   public:
    RangeXYZ<TLattice>() {
        kStart = computeK<TLattice>(0, 0, 0);
        kEnd = computeK<TLattice>(TLattice::subArray[0] - 1, TLattice::subArray[1] - 1, TLattice::subArray[2] - 1) + 1;
    }

    using Iterator = RangeXYZIterator<TLattice, std::array<int, 3>>;
    Iterator begin() { return Iterator(kStart); }
    Iterator end() { return Iterator(kEnd); }

   private:
    int kStart, kEnd;
};

/**
 * \brief Iterator over the x,y,z coordinates and the index k of the local lattice.
 */
template <class TLattice>
class RangeXYZK {
   public:
    RangeXYZK<TLattice>() {
        kStart = computeK<TLattice>(0, 0, 0);
        kEnd = computeK<TLattice>(TLattice::subArray[0] - 1, TLattice::subArray[1] - 1, TLattice::subArray[2] - 1) + 1;
    }

    using Iterator = RangeXYZIterator<TLattice, std::array<int, 4>>;
    Iterator begin() { return Iterator(kStart); }
    Iterator end() { return Iterator(kEnd); }

   private:
    int kStart, kEnd;
};

/**
 * \brief Returns true if the current TLattice point lies on a periodic boundary.
 * \param k Index of current TLattice point.
 * \return True if TLattice point lies on a periodic boundary.
 */
template <class TLattice>
inline bool isPeriodic(int k) {
    int yAtCurrentk = computeY(TLattice::LYdiv, TLattice::LZdiv, k);
    int zAtCurrentk = computeZ(TLattice::LYdiv, TLattice::LZdiv, k);
    int xAtCurrentk = computeX(TLattice::LYdiv, TLattice::LZdiv, k);

    if (TLattice::LZdiv <= 1 || TLattice::LYdiv <= 1 || TLattice::LXdiv <= 1)
        return true;  // If simulation is 2D
    else if (zAtCurrentk == 0 || zAtCurrentk == TLattice::LZdiv - 1)
        return true;  // Edges in Z direction

    else if (yAtCurrentk == 0 || yAtCurrentk == TLattice::LYdiv - 1)
        return true;  // Edges in Y direction

    else if (xAtCurrentk == 0 || xAtCurrentk == TLattice::LXdiv - 1)
        return true;  // Edges in X direction

    return false;
}

template <class TForce>
typename TForce::Method getMethod(TForce& f) {
    return std::declval<typename TForce::Method>();
}

template <class TForce>
typename TForce::Prefactor getForcePrefactor(TForce& f) {
    return std::declval<typename TForce::Prefactor>();
}

void print() {
    if (mpi.rank != 0) return;
    std::cout << std::endl;
}

template <typename T, typename... TArgs>
void print(std::vector<T> first, TArgs... args) {
    if (mpi.rank != 0) return;
    for (auto elem : first) {
        std::cout << elem << " ";
    }
    print(args...);
}

template <typename T, typename... TArgs>
void print(T first, TArgs... args) {
    if (mpi.rank != 0) return;
    std::cout << first << " ";
    print(args...);
}

void printAll() { std::cout << std::endl; }

template <typename T, typename... TArgs>
void printAll(std::vector<T> first, TArgs... args) {
    for (auto elem : first) {
        std::cout << elem << " ";
    }
    printAll(args...);
}

template <typename T, typename... TArgs>
void printAll(T first, TArgs... args) {
    std::cout << first << " ";
    printAll(args...);
}

template <typename T>
std::vector<T> loadTxt(std::string filename, int n = 0) {
    std::vector<T> output;
    if (n != 0) output.reserve(n);

    std::ifstream file(filename);
    if (file.is_open()) {
        T value;
        while (file >> value) {
            output.push_back(value);
        }
        file.close();
    } else {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    return output;
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        return h1 ^ h2;
    }
};

// Compile-time for loop to call func<i>() for i<N
template <int N, typename F>
void constexpr_for(F func) {
    constexpr_for<N>(func, std::make_integer_sequence<int, N>{});
}

template <int N, typename F, int... I>
void constexpr_for(F func, std::integer_sequence<int, I...>) {
    (func.template operator()<I>(), ...);  // Calls func<I>() for each I
}
