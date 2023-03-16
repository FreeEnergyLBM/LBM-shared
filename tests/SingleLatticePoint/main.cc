#include "../../src/Algorithm.hh"
#include "../../src/LBModels.hh"
#include "../../src/Forces.hh"
#include "../../src/Data.hh"
#include "../../src/Stencil.hh"
#include "../../src/Global.hh"
#include "../../src/Service.hh"
#include "../../src/Saving.hh"
#include <iostream>

using SingleComponentWithBodyForce=SingleComponent<D2Q9,Data1,BodyForce,BodyForce>;

using SingleComponentWithoutBodyForce=SingleComponent<D2Q9,Data1>;

using BinaryTest=Binary<D2Q9,Data1>;

int main(){
    data_dir="data/";
    BodyForce Force;
    SingleComponentWithBodyForce test(Force,Force);//get rid of this

    SingleComponentWithoutBodyForce test2;
    BinaryTest test22;

    Algorithm<SingleComponentWithoutBodyForce,BinaryTest> LBM(test2,test22);

    LBM.initialise();
    for (int timestep=0;timestep<=TIMESTEPS;timestep++){
        std::cout<<std::endl<<"############ t="<<timestep<<" ############"<<std::endl;
        LBMPrint(test2);
        LBMPrint(test22);

        LBM.evolve();
        //save(timestep);

    }
    std::cout<<std::endl;
    
    return 0;
}