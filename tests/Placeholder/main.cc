#include "mainInflow.hh"
#include <cstdlib>

int main(int argc, char **argv){

    int seed;
    std::string seedstr;
    mpi.init();

    MPI_Type_create_resized(MPI_INT, 0L, sizeof(Boundary<Lattice::NDIM>), &mMPIBoundary);
    MPI_Type_commit(&mMPIBoundary);

    int blocklengths[3] = {1,1,Lattice::NDIM};
    MPI_Datatype types[3] = {MPI_INT,MPI_CXX_BOOL,MPI_INT8_T};
    std::array<int8_t,Lattice::NDIM> a;
    MPI_Aint offsets[3];
    offsets[0] = offsetof(Boundary<Lattice::NDIM>, Id);
    offsets[1] = offsetof(Boundary<Lattice::NDIM>, IsCorner);
    offsets[2] = offsetof(Boundary<Lattice::NDIM>, NormalDirection) + (size_t)((((char *)(&a)-(char *)(&a[0]))));// + offsetof(std::array<int8_t,TLattice::NDIM>, NormalDirection);
    
    MPI_Type_create_struct(3, blocklengths, offsets, types, &mMPIBoundary);
    MPI_Type_commit(&mMPIBoundary);

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
    //auto flowfield = initFlowField<>();
    auto humidity = initHumidity<>();

    Geometry<Lattice>::initialiseBoundaries(initBoundary,{0,5,6});
    OrderParameter<>::set<Lattice>(initFluid);
    OrderParameterOld<>::set<Lattice>(initFluid);
    //Velocity<>::set<Lattice,Lattice::NDIM,0>(initVelocity);
    //Velocity<>::set<Lattice,Lattice::NDIM,1>(initVelocityY);
    //VelocityOld<>::set<Lattice,Lattice::NDIM,0>(initVelocity);
    //VelocityOld<>::set<Lattice,Lattice::NDIM,1>(initVelocityY);
    Humidity<>::set<Lattice>(initHumidity);
    //Algorithm lbm(flowfield);
    Algorithm lbm(binary,pressure,humidity);//
    //Algorithm lbm(pressure);
    //Algorithm lbm(humidity,binary);
    //Algorithm lbm(binary,pressure);//,humidity);

    SaveHandler<Lattice> saver(datadir);
    saver.saveHeader(timesteps, saveInterval);
    
    for (int timestep=0; timestep<=timesteps; timestep++) {
        TIME=timestep;
        // Save the desired parameters, producing a binary file for each.
        if (timestep%saveInterval==0) {
            if(mpi.rank==0)std::cout<<"Saving at timestep "<<timestep<<"."<<std::endl;
            
            saver.saveBoundaries(timestep);
            saver.saveParameter<Humidity<>>(timestep);
            //saver.SaveParameter<ChemicalPotential<>>(timestep);
            saver.saveParameter<Density<>>(timestep);
            saver.saveParameter<Pressure<>>(timestep);
            saver.saveParameter<OrderParameter<>>(timestep);
            //saver.SaveParameter<LaplacianOrderParameter<>>(timestep);
            saver.saveParameter<MassSink<>>(timestep);
            saver.saveParameter<Velocity<>,Lattice::NDIM>(timestep);
            //saver.SaveParameter<VelocityOld<>,Lattice::NDIM>(timestep);
            //saver.SaveParameter<GradientHumidity<>,Lattice::NDIM>(timestep);
            
        }
        AfterEquilibration(timestep,binary);
        // Evolve by one timestep
        lbm.evolve();
        
    }

}
