#pragma once
#include <vector>
#include <string>
#include <iomanip>
#include <any>
#include <fstream>
#include <iostream>
#include <cstring>
#include "Mpi.hh"
#include "Service.hh"
#include "Parameters.hh"

//Saving.hh: This file contains the functions for saving files.


template<class TLattice>
class SaveHandler {
    public:
        SaveHandler(std::string datadir);
        SaveHandler(SaveHandler& other);

        template<typename... TParameter> void saveVTK(int timestep, TParameter&... params);
        template<typename... TParameter> void saveDAT(int timestep, TParameter&... params);

        inline void saveHeader(int timestep, int saveinterval);

        void saveBoundaries(int timestep);
        void saveBoundariesVTK(int timestep);
        void saveBoundariesDAT(int timestep);

        template<class TParameter, int TNumDir=1> void saveParameter(std::string filename, int instance=-1);
        template<class TParameter, int TNumDir=1> void saveParameter(int timestep, int instance=-1);
        template<class... TParameter> void saveParameter(int timestep, TParameter&... params);

        template<class TParameter, int TNumDir=1> void loadParameter(std::string filename, std::string filetype="bin", int instance=-1);
        template<class TParameter, int TNumDir=1> void loadParameter(int timestep, std::string filetype="bin", int instance=-1);
        template<class... TParameter> void loadParameter(int timestep, TParameter&... params);

        void maskSolid(bool status=true);
        void maskSolid(double value);

    private:
        static bool mMaskSolid;
        static double mMaskValue;
        std::string mDataDir;
        #ifdef MPIPARALLEL
        MPI_Datatype mMPIBoundary;
        #endif
};

template<class TLattice> bool SaveHandler<TLattice>::mMaskSolid = false;
template<class TLattice> double SaveHandler<TLattice>::mMaskValue = 0;


template<class TLattice>
SaveHandler<TLattice>::SaveHandler(std::string datadir) : mDataDir(datadir) {
    int status = system(((std::string)"mkdir -p " + mDataDir).c_str());
    if (status) print("Error creating output directory");
#ifdef MPIPARALLEL
    MPI_Type_create_resized(MPI_INT, 0L, sizeof(Boundary<TLattice::NDIM>), &mMPIBoundary);
    MPI_Type_commit(&mMPIBoundary);
#endif
}

template<class TLattice>
SaveHandler<TLattice>::SaveHandler(SaveHandler& other) : mDataDir(other.datadir) {
    int status = system(((std::string)"mkdir -p " + mDataDir).c_str());
    if (status) print("Error creating output directory");
#ifdef MPIPARALLEL
    MPI_Type_create_resized(MPI_INT, 0L, sizeof(Boundary<TLattice::NDIM>), &mMPIBoundary);
    MPI_Type_commit(&mMPIBoundary);
#endif
}


template<class TLattice>
bool isMasked(int k) {
    int nodeType = Geometry<TLattice>::getBoundaryType(k);
    return (nodeType==1 || nodeType==-1);
}

template<class TLattice, typename T>
void applyMask(std::vector<T>& data, int directions, double maskValue) {
    int nLattice = data.size() / directions;
    for (int k=0; k<nLattice; k++) {
        if (isMasked<TLattice>(k)) {
            for (int iDir=0; iDir<directions; iDir++) {
                int iData = k*directions + iDir;
                data[iData] = maskValue;
            }
        }
    }
}


#ifdef MPIPARALLEL
void writeText(MPI_File file, std::string text) {
    if (mpi.rank == 0) {
        MPI_File_write_shared(file, text.c_str(), text.length(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
}
#else
void writeText(std::ofstream& file, std::string text) {
    file << text;
}
#endif


#ifdef MPIPARALLEL
template<class TLattice, typename T>
void readWriteArray(const char rw, MPI_File file, std::vector<T>& data, int nDim=1, int nInst=1, int instance=-1) {
    int nInstFile = 1;
    if (instance == -1) { // Read all instances
        nInstFile = nInst;
        instance = 0;
    }
    int lsize[5] = {TLattice::subArray[0], TLattice::subArray[1], TLattice::subArray[2], nInstFile, nDim}; // Local array size (to read/write)
    int fsize[5] = {TLattice::LX, TLattice::LY, TLattice::LZ, nInstFile, nDim}; // File array size (global)
    int dsize[5] = {TLattice::LXdiv, TLattice::LYdiv, TLattice::LZdiv, nInst, nDim}; // Data array size (local+halo)
    int flstart[5] = {TLattice::LXMPIOffset, TLattice::LYMPIOffset, TLattice::LZMPIOffset, 0, 0}; // Local array start in file
    int dlstart[5] = {TLattice::HaloXWidth, TLattice::HaloYWidth, TLattice::HaloYWidth, instance, 0}; // Local array start in data

    int nf = fsize[0] * fsize[1] * fsize[2] * fsize[3] * fsize[4];

    // Create subarray objects for writing to file and reading from the data array (or vice versa)
    MPI_Datatype fileSubArray;
    MPI_Datatype dataSubArray;
    MPI_Type_create_subarray(5, fsize, lsize, flstart, MPI_ORDER_C, MPI_DOUBLE, &fileSubArray);
    MPI_Type_create_subarray(5, dsize, lsize, dlstart, MPI_ORDER_C, MPI_DOUBLE, &dataSubArray);
    MPI_Type_commit(&fileSubArray);
    MPI_Type_commit(&dataSubArray);

    // Read / write
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Offset offset;
    MPI_File_get_position_shared(file, &offset);
    MPI_File_set_view(file, offset, MPI_DOUBLE, fileSubArray, "native", MPI_INFO_NULL);
    if (rw == 'r') {
        MPI_File_read_all(file, &data[0], 1, dataSubArray, MPI_STATUS_IGNORE);
    } else if (rw =='w') {
        MPI_File_write_all(file, &data[0], 1, dataSubArray, MPI_STATUS_IGNORE);
    }

    // Move file pointer
    offset += nf * sizeof(T);
    MPI_File_set_view(file, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

    // Delete subarray mpi types
    MPI_Type_free(&fileSubArray);
    MPI_Type_free(&dataSubArray);
}

#else
template<class TLattice, typename T>
void readWriteArray(const char rw, std::fstream& file, std::vector<T>& data, int nDim=1, int nInst=1, int instance=-1) {
    int dInst = (instance==-1) ? 1 : nInst;
    instance = (instance==-1) ? 0 : instance;
    for (int k=0; k<TLattice::N; k++) {
        for (int iInst=instance; iInst<nInst; iInst+=dInst) {
            for (int iDim=0; iDim<nDim; iDim++) {
                int iData = k*nInst*nDim + iInst*nDim + iDim;
                if (rw == 'r') {
                    file.read((char *)(&data[iData]), sizeof(T));
                } else if (rw == 'w') {
                    file.write((char *)(&data[iData]), sizeof(T));
                }
            }
        }
    }
}
#endif


template<typename T>
std::string formatPoint(T data[], std::string fmt, std::string delim1, std::string delim2, int ndim) {
    std::string output = "";
    for (int i=0; i<ndim; i++) {
        char buff[100];
        if constexpr (std::is_same<T, std::string>::value) {
            snprintf(buff, sizeof(buff), fmt.c_str(), data[i].c_str());
        } else {
            snprintf(buff, sizeof(buff), fmt.c_str(), data[i]);
        }
        output += std::string(buff);
        if (i<ndim) output += delim1;
    }
    output += delim2;
    return output;
}


#ifdef MPIPARALLEL
template<class TLattice, typename T>
void writeArrayTxt(MPI_File file, const std::vector<T>& data, std::string fmt, std::string delim="\n", std::string end="NONE", int ndim=1, std::string dimDelim=" ") {// TODO: instances
    int lsize[3] = {TLattice::subArray[0], TLattice::subArray[1], TLattice::subArray[2]};
    int fsize[3] = {TLattice::LX, TLattice::LY, TLattice::LZ};
    int dsize[3] = {TLattice::LXdiv, TLattice::LYdiv, TLattice::LZdiv};
    int flstart[3] = {TLattice::LXMPIOffset, TLattice::LYMPIOffset, TLattice::LZMPIOffset};
    int dlstart[3] = {TLattice::HaloXWidth, TLattice::HaloYWidth, TLattice::HaloYWidth};

    int nf = fsize[0] * fsize[1] * fsize[2];
    int nl = lsize[0] * lsize[1] * lsize[2];

    // Convert data to plaintext
    int numStringLen = formatPoint(&data[0], fmt, dimDelim, delim, ndim).length();
    std::string dataString = "";
    for (int i=dlstart[0]; i<dlstart[0]+lsize[0]; i++) {
        for (int j=dlstart[1]; j<dlstart[1]+lsize[1]; j++) {
            for (int k=dlstart[2]; k<dlstart[2]+lsize[2]; k++) {
                int iData = i*dsize[1]*dsize[2] + j*dsize[2] + k;
                dataString += formatPoint(&data[iData*ndim], fmt, dimDelim, delim, ndim);
            }
        }
    }

    // Create an MPI datatype for writing the text to file
    MPI_Datatype numString;
    MPI_Type_contiguous(numStringLen, MPI_CHAR, &numString);
    MPI_Type_commit(&numString);

    // Create subarray
    MPI_Datatype subarray;
    MPI_Type_create_subarray(3, fsize, lsize, flstart, MPI_ORDER_C, numString, &subarray);
    MPI_Type_commit(&subarray);

    // Get the current offset from the beginning of the file
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Offset offset;
    MPI_File_get_position_shared(file, &offset);

    // Write
    MPI_File_set_view(file, offset, MPI_CHAR, subarray, "native", MPI_INFO_NULL);
    MPI_File_write_all(file, dataString.c_str(), nl, numString, MPI_STATUS_IGNORE);
    MPI_File_set_view(file, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

    // Write any ending and progress the file pointer
    offset += nf * numStringLen;
    if (end != "NONE") {
        if (mpi.rank==0) MPI_File_write_at(file, offset-delim.length(), end.c_str(), end.length(), MPI_CHAR, MPI_STATUS_IGNORE);
        offset += end.length() - delim.length();
    }
    MPI_File_seek_shared(file, offset, MPI_SEEK_CUR);

    // Delete subarray mpi type
    MPI_Type_free(&subarray);
}

#else
template<class TLattice, typename T>
void writeArrayTxt(std::ofstream& file, const std::vector<T>& data, std::string fmt, std::string delim="\n", std::string end="NONE", int ndim=1, std::string dimDelim=" ") {
    // Convert data to plaintext
    std::string dataString = "";
    for (int iData=0; iData<TLattice::N; iData++) {
        dataString += formatPoint(&data[iData*ndim], fmt, dimDelim, delim, ndim);
    }
    file << dataString;
}
#endif

/// MPIPARALLEL part should be corrected!
// #ifdef MPIPARALLEL
// template<class TLattice, typename T>
// void rearrangeArrayTxt(const std::vector<T>& data, std::vector<T>& rearrangedData, int dir) {
//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     // Allocate temporary buffers for each process
//     std::vector<T> localData(data.size() / size);
//     std::vector<T> localRearrangedData(rearrangedData.size() / size);

//     // Distribute data across processes
//     MPI_Scatter(data.data(), data.size() / size, MPI_DOUBLE, localData.data(), data.size() / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//     // Perform local rearrangement
//     for (int i = 0; i < TLattice::LX / size; i++) {
//         for (int j = 0; j < TLattice::LY; j++) {
//             for (int k = 0; k < TLattice::LZ; k++) {
//                 int index = (i * TLattice::LY * TLattice::LZ + j * TLattice::LZ + k) * dir;
//                 int newIndex = (k * TLattice::LY * (TLattice::LX / size) + j * (TLattice::LX / size) + i) * dir;
//                 localRearrangedData[newIndex] = localData[index];
//                 localRearrangedData[newIndex + 1] = localData[index + 1];
//                 localRearrangedData[newIndex + 2] = localData[index + 2];
//             }
//         }
//     }

//     // Gather rearranged data from all processes
//     MPI_Gather(localRearrangedData.data(), localRearrangedData.size(), MPI_DOUBLE, rearrangedData.data(), localRearrangedData.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
// }
// #else
template<class TLattice, typename T>
void rearrangeArrayTxt(const std::vector<T>& data, std::vector<T>& rearrangedData, int dir){

    for (int i=0; i<TLattice::LX; i++) {
            for (int j=0; j<TLattice::LY; j++) { 
                    for (int k=0; k<TLattice::LZ; k++) {    
                int index = (i*TLattice::LY*TLattice::LZ + j*TLattice::LZ + k) * dir;
                int newIndex = (k*TLattice::LY*TLattice::LX + j*TLattice::LX + i) * dir;
                rearrangedData[newIndex] = data[index];
                if (dir != 1)
                {
                    rearrangedData[newIndex + 1] = data[index + 1];
                    rearrangedData[newIndex + 2] = data[index + 2];
                }
            }
        }
    }
}
// #endif

template<class TLattice>
template<typename... TParameter>
void SaveHandler<TLattice>::saveVTK(int timestep, TParameter&... params) {

    std::string filePrefix = "data";
    char filename[512];
    sprintf(filename, "%s/%s_%d.vtk", mDataDir.c_str(), filePrefix.c_str(), timestep);
    #ifdef MPIPARALLEL
        MPI_File file;
        MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    #else
        std::ofstream file(filename);
    #endif

    // Write the header
    char header[1024];
    int nPoints = TLattice::LX * TLattice::LY * TLattice::LZ;
    sprintf(header, "# vtk DataFile Version 3.0\nSimulation data at timestep %d\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA %d\n", timestep, TLattice::LX, TLattice::LY, TLattice::LZ, nPoints);
    writeText(file, header);

    // Write each parameter
    ([&]
    {
        int directions = params.mDirections;
        auto data = params.getParameter();
        // Write header
        char dataHeader[1024];
        if (directions == 1) {
            sprintf(dataHeader, "\nSCALARS %s float\nLOOKUP_TABLE default\n", TParameter::mName);
        } else {
            sprintf(dataHeader, "\nVECTORS %s float\n", TParameter::mName);
        }
        writeText(file, dataHeader);
        // Apply the solid mask, if using
        if (mMaskSolid) applyMask<TLattice>(data, directions, mMaskValue);
        // If 2D vectors, fill the z-component with zeros
        if (directions == 2) {
            directions = 3;
            decltype(data) newData(data.size()/2*3);
            for (int i=0; i<(int)data.size()/2; i++) {
                newData[3*i] = data[2*i];
                newData[3*i+1] = data[2*i+1];
            }
            data = newData;
        }

        // data iterates in the order of z, y, and x, which
        // is not suitable for the visualisation with ParaView.
        // Rearrange the data to reiterate in the order of x, y, and z.
        std::vector<double> rearrangedData(data.size());
        rearrangeArrayTxt<TLattice>(data, rearrangedData, directions);
        
        // Write data
        writeArrayTxt<TLattice>(file, rearrangedData, "%13.8f", "\n", "NONE", directions, " ");
    } (), ...);

    #ifdef MPIPARALLEL
        MPI_File_close(&file);
    #else
        file.close();
    #endif
}


// template<class TLattice>
// template<typename... TParameter>
// void SaveHandler<TLattice>::saveDAT(int timestep, TParameter&... params) {
//     // File layout:
//     // TITLE = [NAME]
//     // VARIABLES = x, y, z, [param]...
//     // ZONE  t= "solid", f= point, I= [LX], J= [LY], K= LZ
//     // [x]\t[y]\t[z]\t[param]...
//     // [x]\t[y]\t[z]\t[param]...
//     // ...

//     std::string filePrefix = "data";
//     char filename[512];
//     sprintf(filename, "%s/%s_%d.dat", mDataDir.c_str(), filePrefix.c_str(), timestep);
//     #ifdef MPIPARALLEL
//         MPI_File file;
//         MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
//     #else
//         std::ofstream file(filename);
//     #endif

//     // Parameter data
//     std::string parameterHeader = "";
//     int nPoints = TLattice::LXdiv * TLattice::LYdiv * TLattice::LZdiv;
//     int nData = 3 * nPoints;
//     ([&] {
//         int directions = params.mDirections;
//         nData += nPoints * directions;
//         if (directions == 1) {
//             parameterHeader += std::string(", ") + TParameter::mName;
//         } else {
//             for (int idir=0; idir<directions; idir++) {
//                 parameterHeader += ", " + (TParameter::mName + std::to_string(idir));
//             }
//         }
//     } (), ...);

//     // Write the header
//     char text[1024];
//     sprintf(text, "TITLE = Simulation data at timestep %d\nVARIABLES = x, y, z%s", timestep, parameterHeader.c_str());
//     writeText(file, text);
//     sprintf(text, "\nZONE  t=\"solid\", f= point, I= %d, J= %d, K= %d\n", TLattice::LX, TLattice::LY, TLattice::LZ);
//     writeText(file, text);

//     // Write the data
//     std::vector<std::string> data;
//     data.reserve(nData);
//     int nWidth = std::to_string(std::max({TLattice::LX, TLattice::LY, TLattice::LZ})).length();
//     for (int x=0; x<TLattice::subArray[0]; x++) {
//         for (int y=0; y<TLattice::subArray[1]; y++) {
//             for (int z=0; z<TLattice::subArray[2]; z++) {
//                 int k = computeK<TLattice>(x, y, z);
//                 std::stringstream pointData;
//                 // Global coordinates
//                 pointData << std::setw(nWidth) << x+TLattice::LXMPIOffset-TLattice::HaloXWidth << "\t";
//                 pointData << std::setw(nWidth) << y+TLattice::LYMPIOffset-TLattice::HaloYWidth << "\t";
//                 pointData << std::setw(nWidth) << z+TLattice::LZMPIOffset-TLattice::HaloZWidth;
//                 // Parameters
//                 ([&] {
//                     int directions = params.mDirections;
//                     for (int idir=0; idir<directions; idir++) {
//                         auto value = params.getParameter()[k*directions+idir];
//                         if (mMaskSolid && isMasked<TLattice>(k)) value = mMaskValue; // Apply the solid mask, if using
//                         pointData << "\t" << std::setw(12) << value;
//                     }
//                 } (), ...);
//                 data.push_back(pointData.str());
//             }
//         }
//     }
//     writeArrayTxt<TLattice>(file, data, "%s");

//     #ifdef MPIPARALLEL
//         MPI_File_close(&file);
//     #else
//         file.close();
//     #endif
// }


template<class TLattice>
template<typename... TParameter>
void SaveHandler<TLattice>::saveDAT(int timestep, TParameter&... params) {
    
    std::string filePrefix = "data";
    char filename[512];
    sprintf(filename, "%s/%s_%d.dat", mDataDir.c_str(), filePrefix.c_str(), timestep);
    #ifdef MPIPARALLEL
        MPI_File file;
        MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    #else
        std::ofstream file(filename);
    #endif

    // Write the header
    char header[1024];
    sprintf(header, "TITLE = \"Simulation data at timestep %d\"\n", timestep);
    writeText(file, header);

    // Write variables
    sprintf(header, "VARIABLES = \"x\", \"y\", \"z\"");
    ([&] {
        int directions = params.mDirections;
        if (directions == 1) {
            sprintf(header + strlen(header), ", \"%s\"", TParameter::mName);
        } else {
            for (int idir = 0; idir < directions; idir++) {
                sprintf(header + strlen(header), ", \"%s_%d\"", TParameter::mName, idir);
            }
        }
    } (), ...);
    strcat(header, "\n");
    writeText(file, header);

    // Write zone information
    sprintf(header, "ZONE T = \"solid\", I = %d, J = %d, K = %d, DATAPACKING = POINT, VARLOCATION = ([3]=CELLCENTERED)\n", TLattice::LX, TLattice::LY, TLattice::LZ);
    writeText(file, header);

    // Write the data
    for (int x = 0; x < TLattice::LX; x++) {
        for (int y = 0; y < TLattice::LY; y++) {
            for (int z = 0; z < TLattice::LZ; z++) {
                int k = computeK<TLattice>(x, y, z);
                // Global coordinates
                writeText(file, std::to_string(x).c_str());
                writeText(file, "\t");
                writeText(file, std::to_string(y).c_str());
                writeText(file, "\t");
                writeText(file, std::to_string(z).c_str());
                writeText(file, "\t");
                // Parameters
                ([&] {
                    int directions = params.mDirections;
                    for (int idir = 0; idir < directions; idir++) {
                        auto value = params.getParameter()[k * directions + idir];
                        if (mMaskSolid && isMasked<TLattice>(k)) value = mMaskValue; // Apply the solid mask, if using
                        writeText(file, std::to_string(value).c_str());
                        writeText(file, "\t");
                    }
                } (), ...);
                writeText(file, "\n");
            }
        }
    }

    #ifdef MPIPARALLEL
        MPI_File_close(&file);
    #else
        file.close();
    #endif
}



template<class TLattice>
inline void SaveHandler<TLattice>::saveHeader(int timestep, int saveinterval) { //Function to save parameter stored in this class

    if(mpi.rank==0){
        print("SAVING HEADER");
        char fdump[256];
        sprintf(fdump, "%s/Header.mat", mDataDir.c_str()); //Buffer containing file name and location.

        std::ofstream fs(fdump, std::ios::out | std::ios::binary);

        fs.write((char *)(&TLattice::LX), sizeof(int));
        fs.write((char *)(&TLattice::LY), sizeof(int));
        fs.write((char *)(&TLattice::LZ), sizeof(int));
        fs.write((char *)(&TLattice::NDIM), sizeof(int));
        fs.write((char *)(&timestep), sizeof(int));
        fs.write((char *)(&saveinterval), sizeof(int));

        fs.close();

    }

}


template<class TLattice>
void SaveHandler<TLattice>::saveBoundaries(int timestep){
    char fdump[512];
    sprintf(fdump, "%s/%s_t%i.mat", mDataDir.c_str(), BoundaryLabels<TLattice::NDIM>::mName, timestep); //Buffer containing file name and location.

#ifdef MPIPARALLEL //When MPI is used we need a different approach for saving as all nodes are trying to write to the file

    MPI_File fh;
    MPI_File_open(MPI_COMM_SELF, fdump, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh); //Open the file using mpi in write only mode
    MPI_File_seek(fh, sizeof(int) * mpi.rank * (TLattice::LX * TLattice::LY * TLattice::LZ) / mpi.size, MPI_SEEK_SET); //Skip to a certain location in the file, currently
    MPI_File_write(fh,&BoundaryLabels<TLattice::NDIM>::template get<TLattice>()[TLattice::HaloSize], (TLattice::N - 2 * TLattice::HaloSize), mMPIBoundary, MPI_STATUSES_IGNORE);
    MPI_File_close(&fh);

#else

    std::ofstream fs(fdump, std::ios::out | std::ios::binary);
    fs.seekp(sizeof(int) * mpi.rank * (TLattice::LX * TLattice::LY * TLattice::LZ) / mpi.size);
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {
        fs.write((char *)(&BoundaryLabels<TLattice::NDIM>::template get<TLattice>()[k].Id), sizeof(int));
    };
    fs.close();

#endif

}

template<class TLattice>
void SaveHandler<TLattice>::saveBoundariesVTK(int timestep) {

    std::string filePrefix = BoundaryLabels<TLattice::NDIM>::mName;
    char filename[512];
    sprintf(filename, "%s/%s_%d.vtk", mDataDir.c_str(), filePrefix.c_str(), timestep);
    #ifdef MPIPARALLEL
        MPI_File file;
        MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    #else
        std::ofstream file(filename);
    #endif

    // Write the header
    char header[1024];
    int nPoints = TLattice::LX * TLattice::LY * TLattice::LZ;
    sprintf(header, "# vtk DataFile Version 3.0\nGeometry Labels\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA %d\n", TLattice::LX, TLattice::LY, TLattice::LZ, nPoints);
    writeText(file, header);

    // Write the labels
    std::vector<int> labels(TLattice::N);
    for (int k = TLattice::HaloSize; k < TLattice::N - TLattice::HaloSize; k++) {
        int value = BoundaryLabels<TLattice::NDIM>::template get<TLattice>()[k].Id;
        labels[k] = value;
    };
    
    // Write header
    char dataHeader[1024];
    sprintf(dataHeader, "\nSCALARS %s float\nLOOKUP_TABLE default\n", BoundaryLabels<TLattice::NDIM>::mName);

    writeText(file, dataHeader);
        
    // data iterates in the order of z, y, and x, which
    // is not suitable for the visualisation with ParaView.
    // Rearrange the data to reiterate in the order of x, y, and z.
    std::vector<int> rearrangedLabels(labels.size());
    rearrangeArrayTxt<TLattice>(labels, rearrangedLabels, 1);
        
    // Write data
    writeArrayTxt<TLattice>(file, rearrangedLabels, "%d", "\n", "NONE", 1, " ");


    #ifdef MPIPARALLEL
        MPI_File_close(&file);
    #else
        file.close();
    #endif
}


template<class TLattice>
void SaveHandler<TLattice>::saveBoundariesDAT(int timestep) {
    
    std::string filePrefix = BoundaryLabels<TLattice::NDIM>::mName;
    char filename[512];
    sprintf(filename, "%s/%s_%d.dat", mDataDir.c_str(), filePrefix.c_str(), timestep);
    #ifdef MPIPARALLEL
        MPI_File file;
        MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    #else
        std::ofstream file(filename);
    #endif

    // Write the header and variable name
    char header[1024];
    sprintf(header, "TITLE = \"Geometry Labels\"\n");
    writeText(file, header);

    sprintf(header, "VARIABLES = \"x\", \"y\", \"z\"");
    sprintf(header + strlen(header), ", \"%s\"", BoundaryLabels<TLattice::NDIM>::mName);
    
    strcat(header, "\n");
    writeText(file, header);

    // Write zone information
    sprintf(header, "ZONE T = \"solid\", I = %d, J = %d, K = %d, DATAPACKING = POINT, VARLOCATION = ([3]=CELLCENTERED)\n", TLattice::LX, TLattice::LY, TLattice::LZ);
    writeText(file, header);

    // Write the data
    for (int x = 0; x < TLattice::LX; x++) {
        for (int y = 0; y < TLattice::LY; y++) {
            for (int z = 0; z < TLattice::LZ; z++) {
                int k = computeK<TLattice>(x, y, z);
                // Global coordinates
                writeText(file, std::to_string(x).c_str());
                writeText(file, "\t");
                writeText(file, std::to_string(y).c_str());
                writeText(file, "\t");
                writeText(file, std::to_string(z).c_str());
                writeText(file, "\t");
                
                // Labels
                auto value = BoundaryLabels<TLattice::NDIM>::template get<TLattice>()[k].Id;
                writeText(file, std::to_string(value).c_str());
                writeText(file, "\t");
               
                writeText(file, "\n");
            }
        }
    }

    #ifdef MPIPARALLEL
        MPI_File_close(&file);
    #else
        file.close();
    #endif
}

template<class TLattice>
template<class TParameter, int TNumDir>
void SaveHandler<TLattice>::saveParameter(std::string filename, int instance) {
    // Setup array for saving
    std::vector<typename TParameter::ParamType>& param = TParameter::template get<TLattice,TNumDir>();
    if (mMaskSolid) applyMask<TLattice>(param, TNumDir, mMaskValue);

    // Open and write file
    #ifdef MPIPARALLEL
        MPI_File file;
        MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    #else
        std::fstream file(filename.c_str(), std::fstream::out);
    #endif

    readWriteArray<TLattice>('w', file, param, TNumDir, TParameter::instances, instance);

    #ifdef MPIPARALLEL
        MPI_File_close(&file);
    #else
        file.close();
    #endif
}

template<class TLattice>
template<class TParameter, int TNumDir>
void SaveHandler<TLattice>::saveParameter(int timestep, int instance) {
    char filename[512];
    if (instance == -1) {
        sprintf(filename, "%s/%s_t%d.mat", mDataDir.c_str(), TParameter::mName, timestep);
    } else {
        sprintf(filename, "%s/%s%d_t%d.mat", mDataDir.c_str(), TParameter::mName, instance, timestep);
    }
    saveParameter<TParameter,TNumDir>(filename, instance);
}

template<class TLattice>
template<class... TParameter>
void SaveHandler<TLattice>::saveParameter(int timestep, TParameter&... params) {
    (saveParameter<TParameter>(timestep),...);
}



template<class TLattice>
void SaveHandler<TLattice>::maskSolid(bool status) {
    mMaskSolid = status;
}


template<class TLattice>
void SaveHandler<TLattice>::maskSolid(double value) {
    mMaskSolid = true;
    mMaskValue = value;
}


//==== Loading ====//
template<class TLattice>
template<class TParameter, int TNumDir>
void SaveHandler<TLattice>::loadParameter(std::string filename, std::string filetype, int instance) {
    std::vector<typename TParameter::ParamType> &param = TParameter::template get<TLattice,TNumDir>();

    if (filetype == "bin") {
        // Open and read file
        #ifdef MPIPARALLEL
            MPI_File file;
            MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
        #else
            std::fstream file(filename.c_str());
        #endif
        readWriteArray<TLattice>('r', file, param, TNumDir, TParameter::instances, instance);
        #ifdef MPIPARALLEL
            MPI_File_close(&file);
        #else
            file.close();
        #endif

    } else if (filetype == "txt") {
        std::ifstream fs(filename.c_str());
        int xStart = TLattice::LXMPIOffset - TLattice::HaloXWidth;
        int yStart = TLattice::LYMPIOffset - TLattice::HaloYWidth;
        int zStart = TLattice::LZMPIOffset - TLattice::HaloZWidth;
        int xStop = std::min(xStart + TLattice::LXdiv, TLattice::LX);
        int yStop = std::min(yStart + TLattice::LYdiv, TLattice::LY);
        int zStop = std::min(zStart + TLattice::LZdiv, TLattice::LZ);
        int dInst = (instance==-1) ? 1 : TParameter::instances;
        if (instance==-1) instance = 0;
        double dummy;
        for (int xg=0; xg<xStop; xg++) {
            for (int yg=0; yg<yStop; yg++) {
                for (int zg=0; zg<zStop; zg++) {
                    int k = computeKFromGlobal<TLattice>(xg%TLattice::LX, yg%TLattice::LY, zg%TLattice::LZ);
                    for (int iInst=instance; iInst<TParameter::instances; iInst+=dInst) {
                        for (int iDim=0; iDim<TNumDir; iDim++) {
                            if (k == -1) {
                                fs >> dummy;
                            } else {
                                int i = k*TParameter::instances*TNumDir + iInst*TNumDir + iDim;
                                fs >> param[i];
                            }
                        }
                    }
                }
            }
        }
        fs.close();
    }

    // Ensure the parameter is set to initialised
    std::map<int,bool> &initialised = TParameter::template getInstance<TLattice,TNumDir>().mmInitialised;
    for (size_t i=0; i<param.size(); i++) {
        initialised[i] = true;
    }
}


template<class TLattice>
template<class TParameter, int TNumDir>
void SaveHandler<TLattice>::loadParameter(int timestep, std::string filetype, int instance) {
    std::string extension = (filetype=="bin") ? "mat" : filetype;
    char filename[512];
    if (instance==-1) {
        sprintf(filename, "%s/%s_t%i.%s", mDataDir.c_str(), TParameter::mName, timestep, extension.c_str());
    } else {
        sprintf(filename, "%s/%s%d_t%i.%s", mDataDir.c_str(), TParameter::mName, instance, timestep, extension.c_str());
    }
    loadParameter<TParameter,TNumDir>(filename, filetype, instance);
}

template<class TLattice>
template<class... TParameter>
void SaveHandler<TLattice>::loadParameter(int timestep, TParameter&... params) {
    (loadParameter<TParameter>(timestep),...);
}
