#include "main.hh"

int main(int argc, char **argv){

    mpi.init();
    initParams("input.txt");

    auto binary = initBinary<>();
    auto pressure = initPressure<>();
    auto humidity = initHumidity<>();

    Geometry<Lattice>::initialiseBoundaries(initBoundary);
    OrderParameter<>::set<Lattice>(initFluid);
    OrderParameterOld<>::set<Lattice>(initFluid);
    //Humidity<>::set<Lattice>(initHumidity);

    //Algorithm lbm(pressure,binary,humidity);
    Algorithm lbm(binary,pressure);//,humidity);

    ParameterSave<Lattice> saver(datadir);
    saver.SaveHeader(timesteps, saveInterval);
    
    for (int timestep=0; timestep<=timesteps; timestep++) {

        // Save the desired parameters, producing a binary file for each.
        if (timestep%saveInterval==0) {
            if(mpi.rank==0)std::cout<<"Saving at timestep "<<timestep<<"."<<std::endl;
            
            saver.SaveBoundaries(timestep);
            saver.SaveParameter<Humidity<>>(timestep);
            saver.SaveParameter<ChemicalPotential<>>(timestep);
            saver.SaveParameter<Density<>>(timestep);
            saver.SaveParameter<Pressure<>>(timestep);
            saver.SaveParameter<OrderParameter<>>(timestep);
            saver.SaveParameter<LaplacianOrderParameter<>>(timestep);
            saver.SaveParameter<MassSink<>>(timestep);
            saver.SaveParameter<Velocity<>,Lattice::NDIM>(timestep);
            saver.SaveParameter<VelocityOld<>,Lattice::NDIM>(timestep);
            saver.SaveParameter<GradientHumidity<>,Lattice::NDIM>(timestep);
            
        }
        
        // Evolve by one timestep
        lbm.evolve();
        
    }

}