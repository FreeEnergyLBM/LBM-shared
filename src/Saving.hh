#ifndef SAVING_HEADER
#define SAVING_HEADER
#include <fstream>
#include <iostream>
std::string data_dir;
void save(int t){//THIS IS TEMPORARY
    char fdump[512];
    sprintf(fdump, (data_dir+"conf_k_t%li.mat").c_str(),
    t+10000000000
	  );
    Density<double> m_Density;
    Velocity<double,NDIM> m_Velocity;
    std::ofstream fs(fdump, std::ios::out | std::ios::binary );

    for (int k = 0; k < N; k++ ) { 

        fs.write((char *)(&m_Density.getParameter(k)), sizeof(double));
        fs.write((char *)(&m_Velocity.getParameter(k*3)), sizeof(double));
        fs.write((char *)(&m_Velocity.getParameter(k*3+1)), sizeof(double));
        fs.write((char *)(&m_Velocity.getParameter(k*3+2)), sizeof(double));

	  };
}
#endif