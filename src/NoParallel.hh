/**TDistribution&
 * \brief NoParallel is a dummy implementation of Parallel that does not split the lattice or perform any communication.
 */
class NoParallel : public Parallel<NoParallel, 0> {
   public:
    // temporary implementation before a better idea
    static constexpr int mNumDirections = 0;       // number of communication directions in 1D parallelisation
    static constexpr int mCommDirection[1] = {0};  // communication direction along X

    template <class TLattice>
    void init(){};

    template <class TLattice, class TStencil>
    static void createDistributionType() {
#ifdef MPIPARALLEL
        NoParallel::mDistributionType = MPI_DOUBLE;
#endif
    };
};
