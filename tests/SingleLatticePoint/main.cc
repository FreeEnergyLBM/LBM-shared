#include "../../src/Algorithm.hh"
#include "../../src/LBModels.hh"
#include "../../src/Forces.hh"
#include "../../src/Data.hh"
#include "../../src/Stencil.hh"
#include "../../src/Global.hh"
#include <iostream>

using SingleComponentWithBodyForce=SingleComponent<D2Q9,Data_Placeholder,BodyForce>;
using SingleComponentWithoutBodyForce=SingleComponent<D2Q9,Data_Placeholder>;

int main(){
    BodyForce Force;
    SingleComponentWithBodyForce test(Force);
    SingleComponentWithoutBodyForce test2;
    Algorithm<SingleComponentWithoutBodyForce> LBM(test2);
    LBM.initialise();
    for (int timestep=0;timestep<TIMESTEPS;timestep++){
        LBM.evolve();
        std::cout<<" t="<<timestep<<std::endl;
        for(int i=0;i<9;i++)std::cout<<test2.getDistribution()[i]<<" "<<std::flush;
    }
    
    
    return 0;
}

    