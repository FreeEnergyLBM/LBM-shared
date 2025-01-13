#include <math.h>

#include <lbm.hh>

// This example demonstrates how to couple the binary model with the porous model
// to simulate capillary driven flow from/to a porous medium.

const int lx = 100;  // Size of domain in x direction
const int ly = 5;    // Size of domain in y direction

const int timeSteps = 500;    // Number of time steps to run the simulation
const int saveInterval = 50;  // Interval to save global data

// Parameters to control the surface tension and width of the diffuse interface
// Surface tension (in lattice units) = sqrt(8*kappa*A/9)
// Interface width (in lattice units) = sqrt(kappa/A)
const double binaryA = 1.0e-3;
const double binaryKappa = 2 * binaryA;
const double contactAngle = 90;  // Contact angle of the liquid on the solid

// Relaxation times of each component. tau1 corresponds to phi=1.0, tau2 corresponds to phi=-1.0
// Viscosity (in lattice units) = 1.0/3.0 * (tau - 0.5)
double tau1 = 1.0;
double tau2 = 1.0;

const double porosity = 0.7;         // Porosity of PolyCrystalline Salt
const double permeability = 2286.0;  // Permeability of PolyCrystalline Salt

// Irreducible wetting and non-wetting saturations
double Sw_ir = 0.4;
double Snw_ir = 0.0;

// Capillary pressure parameters
double P_0 = 1.0E-2;
double P_c_inf = 10 * P_0;
double lambda = 2.0;

// Set up the lattice, including the resolution and data/parallelisation method
using Lattice = LatticeProperties<ParallelX<1>, lx, ly>;

// Function used to define the geometry
int initSolid(int k) {
    //-1: BulkSolid = -1, Fluid = 0, Wall = 1, InletWall = 3, OutletWall = 4
    auto [x, y, z] = computeXYZ<Lattice>(k);
    if (x >= lx / 2 - 1 && x <= lx / 2) {
        return 7;
    } else {
        return 0;
    }
}

// Function used to define the fluid
double initFluid(int k) {
    auto [x, y, z] = computeXYZ<Lattice>(k);
    if (initSolid(k) == 7) {
        return 0.0;
    } else {
        if (x < lx / 2 - 1) {
            return 1.0;
        } else {
            return -1.0;
        }
    }
}

// Function used to initialise the porosity
double initPorosity(const int k) { return porosity; }

// Function used to initialise the permeability
double initPermeability(const int k) { return permeability; }

// Function used to initialise the saturation
double initSaturation(const int k) {
    auto [x, y, z] = computeXYZ<Lattice>(k);
    if (x >= lx / 2 - 1 && x <= lx / 2) {
        return 1.0;
    } else {
        return 0.0;
    }
}

int main(int argc, char **argv) {
    mpi.init();

    using TraitFlowFieldBinary = DefaultTraitFlowFieldBinary<Lattice>::SetCollisionOperator<MRT>;

    // MassExchangeCalculator must be used in the Binary model (not porous) to be able to run MassExchangeSource for the
    // binary and porous models
    using TraitBinary = DefaultTraitBinary<Lattice>::template SetProcessor<
        std::tuple<GradientsMultiStencil<OrderParameter<>, CentralXYZ, LaplacianCentral>>,
        std::tuple<ChemicalPotentialCalculatorBinary, CubicWetting>, std::tuple<ResetMassExchangeParameter>,
        std::tuple<MassExchangeCalculator>>::AddForce<MassExchangeSource<EvaporationSourceMethod>>;

    using TraitWettingPorousPoiseuille = DefaultTraitMultiphasePorousFlowFieldInSalt<Lattice>::SetProcessor<
        std::tuple<UpdateRelPermAndPc>,
        std::tuple<GradientsMultiParam<CentralXYZ, Saturation<>, CapillaryPressure<>>>>::
        AddForce<CapillaryForcePorous<GuoMultiphasePorous, Gradient, WettingRelativePermeability<>, VelocityPorous<>>,
                 EvaporationPhaseSource<EvaporationSourceMethod>, MassExchangeSource<EvaporationSourceMethod>>;

    // Define the models to be used
    FlowFieldBinary<Lattice, TraitFlowFieldBinary>
        flowFieldModel;  // Flowfield (navier stokes solver) that can be used with the binary model
    Binary<Lattice, TraitBinary> componentSeparationModel;  // Binary model with hybrid equilibrium and forcing term

    MultiphasePorousFlowFieldInSalt<Lattice, TraitWettingPorousPoiseuille> saltCrystalModel;

    // Pass the relaxation times to each model
    flowFieldModel.setTau1(tau1);
    flowFieldModel.setTau2(tau2);

    componentSeparationModel.setTau1(tau1);
    componentSeparationModel.setTau2(tau2);
    componentSeparationModel.setA(binaryA);

    // Set the relaxation time to the default value
    saltCrystalModel.setTau(1.0);

    // Pass the surface tension/interface width parameters to the relevant preprocessor
    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setA(binaryA);
    componentSeparationModel.getProcessor<ChemicalPotentialCalculatorBinary>().setKappa(binaryKappa);

    // Define the solid boundaries
    Geometry<Lattice>::initialiseBoundaries(initSolid);

    // Set the collide ID for each model
    flowFieldModel.setCollideID({0});
    componentSeparationModel.setCollideID({0});
    saltCrystalModel.setCollideID({7});

    // Set the bounce-back boundary conditions for all models
    flowFieldModel.getBoundary<BounceBack>().setNodeID({1, 7});
    componentSeparationModel.getBoundary<BounceBack>().setNodeID({1, 7});
    saltCrystalModel.getBoundary<BounceBackHLBM>().setNodeID({1, 0});

    componentSeparationModel.getProcessor<CubicWetting>().setNodeID({1, 7});
    componentSeparationModel.getProcessor<CubicWetting>().setThetaDegrees(contactAngle);
    componentSeparationModel.getProcessor<CubicWetting>().setAlpha(sqrt(2));
    componentSeparationModel.getProcessor<CubicWetting>().useSinglePhaseCheck = true;

    // Set the irreducible wetting and non-wetting saturations
    saltCrystalModel.template getProcessor<UpdateRelPermAndPc>().setSwIr(Sw_ir);
    saltCrystalModel.template getProcessor<UpdateRelPermAndPc>().setSnwIr(Snw_ir);

    // Set the capillary pressure parameters
    saltCrystalModel.template getProcessor<UpdateRelPermAndPc>().setP0(P_0);
    saltCrystalModel.template getProcessor<UpdateRelPermAndPc>().setPcInf(P_c_inf);
    saltCrystalModel.template getProcessor<UpdateRelPermAndPc>().setLambda(lambda);

    // Set the Leverett J-function parameters
    saltCrystalModel.template getProcessor<UpdateRelPermAndPc>().setK0(permeability);
    saltCrystalModel.template getProcessor<UpdateRelPermAndPc>().setPhi0(porosity);

    // Set the relative permeability exponents
    saltCrystalModel.template getProcessor<UpdateRelPermAndPc>().setN(1.0);

    // Initialise the fluid using the function above
    OrderParameter<>::set<Lattice>(initFluid);

    // Initialise the porosity using the function above
    Porosity<>::set<Lattice>(initPorosity);

    // Initialise the permeability using the function above
    Permeability<>::set<Lattice>(initPermeability);

    // Initialise the saturation using the function above
    Saturation<>::set<Lattice>(initSaturation);

    // Algorithm creates an object that can run our chosen LBM model
    Algorithm lbm(flowFieldModel, componentSeparationModel, saltCrystalModel);

    // Set up the handler object for saving data
    SaveHandler<Lattice> saver("data/");
    saver.maskSolid();

    // int timestep = 0;
    saver.saveBoundariesVTK(0);

    // Perform the main LBM loop

    for (int timestep = 0; timestep < timeSteps; timestep++) {
        // Save the desired parameters, producing a binary file for each.
        if (timestep % saveInterval == 0) {
            std::cout << "Saving at timestep " << timestep << "." << std::endl;
            saver.saveVTK(timestep, Density<>::template getInstance<Lattice>(),
                          OrderParameter<>::template getInstance<Lattice>(),
                          Saturation<>::template getInstance<Lattice>(),
                          Velocity<>::template getInstance<Lattice, Lattice::NDIM>(),
                          VelocityPorous<>::template getInstance<Lattice, Lattice::NDIM>());
        }

        // Evolve by one timestep
        lbm.evolve();
    }

    std::cout << "Simulation complete." << std::endl;
    return 0;
}
