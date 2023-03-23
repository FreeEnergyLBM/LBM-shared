#include "../../src/Algorithm.hh"
#include "../../src/LBModels/Models.hh"
#include "../../src/Forces/Forces.hh"
#include "../../src/BoundaryModels/Boundaries.hh"
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

//TODO BEFORE HACKATHON
//BASIC MPI
//BINARY
//BOUNCEBACK
//COMMENTS
//CLEANUP (SAVING, EXCEPTIONS etc)
//CLEANUP OLD CODE

struct traitF{
    using Stencil=D2Q9;
    template<class model>
    using Data=Data1<Stencil,model>;
    using Boundaries=std::tuple<BounceBack&>;
    using Forces=std::tuple<BodyForce&>;
};
struct traitG{
    using Stencil=D2Q9;  
    template<class model>
    using Data=Data1<Stencil,model>;
    using Boundaries=std::tuple<BounceBack&>;
    using Forces=std::tuple<>;
};

int main(){
    data_dir="data/";
    
    BodyForce Force;
    auto ForcesF=GenerateTuple(Force);
    BounceBack Boundary;
    auto BoundaryF=GenerateTuple(Boundary);
    //SingleComponent<traits> test(Force,Force);//get rid of this

    SingleComponent<traitF> test2(ForcesF,BoundaryF);
    Binary<traitG> test22(BoundaryF);

    Algorithm<SingleComponent<traitF>,Binary<traitG>> LBM(test2,test22);

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