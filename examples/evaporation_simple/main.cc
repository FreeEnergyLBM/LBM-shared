#include <lbm.hh>

// This example performs a simple evaporation where the evaporation rate is proportional to the surface area (using
// SimpleMassLossCalculator)

const int lx = 60;
const int ly = 20;
const int lz = 1;

const int timesteps = 100000;
const int saveInterval = 1000;

const double binaryA = 0.0025;
const double binaryKappa = 0.005;
const double contactAngle = 60;
const double x0 = 2 / 3. * lx;

using Lattice = LatticeProperties<ParallelX<1>, lx, ly, lz>;

int initSolid(const int k) {
    int y = computeY<Lattice>(k);
    int x = computeX<Lattice>(k);
    if (x < 2 || x >= lx - 2) return 1;
    if (y < 2 || y >= ly - 2)
        return 1;
    else
        return 0;
}

double initFluid(int k) {
    int x = computeX<Lattice>(k);
    double width = sqrt(8 * binaryKappa / binaryA);
    return 0.5 * (1 - tanh(2 * (x - x0) / width));
}

using OrderParameterGradients =
    GradientsMultiStencil<OrderParameter<>, CentralXYZ, CentralQ, MixedXYZ, MixedQ, LaplacianCentralWetting>;
using ChemicalPotentialGradients = Gradients<ChemicalPotential<>, LaplacianCentralMirrorSolid>;
using TraitBinary =
    DefaultTraitBinaryLeeHumidity<Lattice>::SetProcessor<ChemicalPotentialGradients, OrderParameterGradients,
                                                         ChemicalPotentialCalculatorBinaryLee, SimpleMassLossCalculator,
                                                         NoFluxSolid<OrderParameter<>>>;

int main(int argc, char **argv) {
    mpi.init();

    // Define the model
    BinaryLeeHumidity<Lattice, TraitBinary> binary;
    binary.getProcessor<ChemicalPotentialCalculatorBinaryLee>().setA(binaryA);
    binary.getProcessor<ChemicalPotentialCalculatorBinaryLee>().setKappa(binaryKappa);

    double wettingPrefactor = -cos(contactAngle * M_PI / 180.0) * sqrt(2 * binaryA / binaryKappa);
    binary.getProcessor<OrderParameterGradients>().setWettingPrefactor(wettingPrefactor);
    binary.getProcessor<OrderParameterGradients>().setBoundaryID(1);
    binary.getProcessor<ChemicalPotentialGradients>().setBoundaryID(1);

    // Boundaries
    Geometry<Lattice>::initialiseBoundaries(initSolid);
    binary.getBoundary<BounceBack>().setNodeID(1);
    binary.getProcessor<NoFluxSolid<OrderParameter<>>>().setNodeID(1);

    // Initialise
    OrderParameter<>::set<Lattice>(initFluid);

    SaveHandler<Lattice> saver("data/");
    saver.saveHeader(timesteps, saveInterval);
    saver.maskSolid();

    // Main loop
    Algorithm lbm(binary);
    for (int timestep = 0; timestep <= timesteps; timestep++) {
        if (timestep % saveInterval == 0) {
            print("Saving at timestep:", timestep);
            // saver.saveDAT(timestep, OrderParameter<>::template getInstance<Lattice>());
            saver.saveVTK(timestep, OrderParameter<>::template getInstance<Lattice>(),
                          MassSink<>::template getInstance<Lattice>());
        }
        lbm.evolve();
    }

    return 0;
}
