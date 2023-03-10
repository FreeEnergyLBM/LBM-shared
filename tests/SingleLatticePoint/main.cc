#include "../../src/Algorithm.hh"
#include "../../src/LBModels.hh"
#include "../../src/Forces.hh"
#include "../../src/Data.hh"
#include "../../src/Stencil.hh"
#include "../../src/Global.hh"
#include "../../src/Service.hh"
#include <iostream>

using SingleComponentWithBodyForce=SingleComponent<D3Q19,Data1,BodyForce>;

using SingleComponentWithoutBodyForce=SingleComponent<D3Q19,Data1>;

int main(){

    BodyForce Force;
    SingleComponentWithBodyForce test(Force);//get rid of this

    SingleComponentWithoutBodyForce test2;

    Algorithm<SingleComponentWithoutBodyForce> LBM(test2);
    
    LBM.initialise();
    std::cout<<std::endl<<"############ t="<<-1<<" ############"<<std::endl;
    LBMPrint(test2);
    for (int timestep=0;timestep<TIMESTEPS;timestep++){
        LBM.evolve();
        std::cout<<std::endl<<"############ t="<<timestep<<" ############"<<std::endl;
        LBMPrint(test2);
    }
    std::cout<<std::endl;
    
    return 0;
}