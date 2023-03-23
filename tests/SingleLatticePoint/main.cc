#include "../../src/Algorithm.hh"
#include "../../src/LBModels/Models.hh"
#include "../../src/Forces/Forces.hh"
#include "../../src/Data.hh"
#include "../../src/Stencil.hh"
#include "../../src/Global.hh"
#include "../../src/Service.hh"
#include "../../src/Saving.hh"
#include <iostream>

//Code wishlist:

//	Grid refinement
//	OpenMP
//	Advanced data - think about affinity
//	DDF shifting
//  Use BLAS
// Ci script
// Debugger
// Sphinx

//USE TRAITS INSTEAD, EG
//struct trait{
// typedef D2Q9 Stencil;    
// typedef Data1 Data;
// typedef std::tuple<BodyForce> Forces;
// ...  
//}
//
//SingleComponent<trait> Modelwithtraits; 

//TODO BEFORE HACKATHON
//BASIC MPI
//BINARY
//BOUNCEBACK
//COMMENTS
//CLEANUP (SAVING, EXCEPTIONS etc)
//CLEANUP OLD CODE

int main(){
    data_dir="data/";
    BodyForce Force;
    //SingleComponent<traits> test(Force,Force);//get rid of this

    SingleComponent<trait<D2Q9,Data1,BodyForce>> test2(Force);
    Binary<trait<D2Q9,Data1>> test22;

    Algorithm<SingleComponent<trait<D2Q9,Data1,BodyForce>>,Binary<trait<D2Q9,Data1>>> LBM(test2,test22);

    LBM.initialise();
    for (int timestep=0;timestep<=TIMESTEPS;timestep++){
        std::cout<<std::endl<<"############ t="<<timestep<<" ############"<<std::endl;
        LBMPrint(test2);
        LBMPrint(test22);

        LBM.evolve();

    }
    std::cout<<std::endl;
    
    return 0;
}