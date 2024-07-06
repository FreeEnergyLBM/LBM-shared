/**
 * \brief ParallelZ contains functions and data for MPI parallelisation divided evenly in the Z direction.
 * This class contains functions to communicate parameters for gradient calculations and to communicate
 * distributions for streaming.
 */
template <int TWidth = 1>
class ParallelZ : public Parallel<ParallelZ<TWidth>, TWidth> {
   public:
    static constexpr int mNumDirections = 1;  // number of communication directions in 1D parallelisation
    static constexpr int mCommDirection[mNumDirections] = {2};  // communication direction along Z

    /**
     * \brief Initialise MPI variables for the parallelisation method.
     */
    template <class TLattice>
    void init();

    template <class TLattice, class TStencil>
    static void createDistributionType();  //!< Create datatype for streaming distributions
};

template <int TWidth>
template <class TLattice>
void ParallelZ<TWidth>::init() {
    if (mpi.size <= 1) return;

    const int faceSize = TLattice::Face[0] = TLattice::LX * TLattice::LY;

    // Split lattice
    TLattice::HaloZWidth = TLattice::Width = this->mMaxWidth;
    TLattice::HaloSize = faceSize * this->mMaxWidth;

    if (TLattice::LZ % mpi.size == 0) {
        TLattice::LZdiv = (TLattice::LZ / mpi.size + 2 * TWidth);
        TLattice::LZMPIOffset = (TLattice::LZdiv - 2 * TWidth) * mpi.rank;
        TLattice::subArray[2] = TLattice::LZ / mpi.size;
    } else {
        std::string err_message =
            "Currently, the size of the domain in the Z direction must be divisible by the number of mpi ranks.";
        throw std::runtime_error(err_message);
    }
    TLattice::N = TLattice::LZdiv * faceSize;

    // Define neighbor processors
    int frontNeighbor = mpi.rank - 1;
    if (frontNeighbor == -1) frontNeighbor = mpi.size - 1;
    int backNeighbor = mpi.rank + 1;
    if (backNeighbor == mpi.size) backNeighbor = 0;
    this->mNeighbors = {frontNeighbor, backNeighbor};

    // Create communication objects
    this->mI0Send = std::vector<int>(2);
    this->mI0Send[0] = 0;
    this->mI0Send[1] = TWidth * faceSize;
    this->mI0Recv = std::vector<int>(2);
    this->mI0Recv[0] = 2 * TWidth * faceSize;
    this->mI0Recv[1] = 3 * TWidth * faceSize;

    this->mI0SendDistr = std::vector<int>(2);
    this->mI0SendDistr[0] = 0;
    this->mI0SendDistr[1] = faceSize;  // 3*faceSize
    this->mI0RecvDistr = std::vector<int>(2);
    this->mI0RecvDistr[0] = 2 * faceSize;  // faceSize
    this->mI0RecvDistr[1] = 3 * faceSize;  // 2*faceSize

// Create MPI buffer
#ifdef MPIPARALLEL
    const int bufSize = (faceSize * TWidth * (5 + 2) * 2 * 2 + 1000) * sizeof(double);

    if (bufSize > MPIBUFFERSIZE) {
        if (MPIBUFFERSIZE != 0) MPI_Buffer_detach(MPIBUFFER, &MPIBUFFERSIZE);
        if (MPIBUFFERSIZE != 0) delete[] MPIBUFFER;
        MPIBUFFER = new char[bufSize];
        MPI_Buffer_attach(MPIBUFFER, bufSize);
        MPIBUFFERSIZE = bufSize;
    }
#endif
}

template <int TWidth>
template <class TLattice, class TStencil>
void ParallelZ<TWidth>::createDistributionType() {
#ifdef MPIPARALLEL
    const int faceSize = TLattice::LX * TLattice::LY;
    MPI_Type_vector(faceSize, 1, TStencil::Q, mpi_get_type<double, TLattice>(), &ParallelZ<TWidth>::mDistributionType);
    MPI_Type_commit(&ParallelZ<TWidth>::mDistributionType);
    ParallelZ<TWidth>::isCustomDistributionTypeCreated = true;
#endif
}
