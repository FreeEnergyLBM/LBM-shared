#include "main.hh"
#include <cstdlib>

int main(int argc, char **argv){

    int seed;
    std::string seedstr;
    mpi.init();
    //exit(1);
    if (argc!=0){
        //std::cerr<<argv[1]<<std::endl;
        seed=std::atoi(argv[1]);
        //std::cerr<<seed<<std::endl;
        seedstr=std::to_string(seed);
        initParams("input/input"+seedstr+".txt");
    }
    else{
        seedstr="";
        initParams("input.txt");
    }

    if(mpi.rank==0&&argc==0){
        int ret;
        std::string tmp="rm input/input"+seedstr+".txt";
        const char *array = tmp.c_str();
        ret = std::system(array);
    }

    auto binary = initBinary<>();
    auto pressure = initPressure<>();
    auto humidity = initHumidity<>();

    Geometry<Lattice>::initialiseBoundaries(initBoundary,{0,5,6});
    OrderParameter<>::set<Lattice>(initFluid);
    OrderParameterOld<>::set<Lattice>(initFluid);
    Humidity<>::set<Lattice>(initHumidity);

    Algorithm lbm(binary,pressure,humidity);//
    //Algorithm lbm(humidity);//,pressure,binary);
    //Algorithm lbm(binary,pressure);//,humidity);

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
        AfterEquilibration(timestep,binary);
        // Evolve by one timestep
        lbm.evolve();
        
    }

}