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

        template<class TParameter, int TNumDir=1> void saveParameter(int timestep);
        template<class... TParameter> void saveParameter(int timestep, TParameter&... params);

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
void writeArray(MPI_File file, std::vector<T> data) {
    int lsize[3] = {TLattice::subArray[0], TLattice::subArray[1], TLattice::subArray[2]};
    int gsize[3] = {TLattice::LX, TLattice::LY, TLattice::LZ};
    int dsize[3] = {TLattice::LXdiv, TLattice::LYdiv, TLattice::LZdiv};
    int glstart[3] = {TLattice::LXMPIOffset, TLattice::LYMPIOffset, TLattice::LZMPIOffset};
    int dlstart[3] = {TLattice::HaloXWidth, TLattice::HaloYWidth, TLattice::HaloYWidth};

    int ng = gsize[0] * gsize[1] * gsize[2];

    // Create subarray objects for writing to file and reading from the data array
    MPI_Datatype fileSubArray;
    MPI_Datatype dataSubArray;
    MPI_Type_create_subarray(2, gsize, lsize, glstart, MPI_ORDER_C, MPI_FLOAT, &fileSubArray);
    MPI_Type_create_subarray(2, dsize, lsize, dlstart, MPI_ORDER_C, MPI_FLOAT, &dataSubArray);
    MPI_Type_commit(&fileSubArray);
    MPI_Type_commit(&dataSubArray);

    // Write
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Offset offset;
    MPI_File_get_position_shared(file, &offset);
    MPI_File_set_view(file, offset, MPI_FLOAT, fileSubArray, "native", MPI_INFO_NULL);
    MPI_File_write_all(file, &data[0], 1, dataSubArray, MPI_STATUS_IGNORE);

    // Move file pointer
    offset += ng * sizeof(data[0]);
    MPI_File_set_view(file, offset, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
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
void writeArrayTxt(MPI_File file, std::vector<T> data, std::string fmt, std::string delim="\n", std::string end="NONE", int ndim=1, std::string dimDelim=" ") {
    int lsize[3] = {TLattice::subArray[0], TLattice::subArray[1], TLattice::subArray[2]};
    int gsize[3] = {TLattice::LX, TLattice::LY, TLattice::LZ};
    int dsize[3] = {TLattice::LXdiv, TLattice::LYdiv, TLattice::LZdiv};
    int glstart[3] = {TLattice::LXMPIOffset, TLattice::LYMPIOffset, TLattice::LZMPIOffset};
    int dlstart[3] = {TLattice::HaloXWidth, TLattice::HaloYWidth, TLattice::HaloYWidth};

    int ng = gsize[0] * gsize[1] * gsize[2];
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
    MPI_Type_create_subarray(2, gsize, lsize, glstart, MPI_ORDER_C, numString, &subarray);
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
    offset += ng * numStringLen;
    if (end != "NONE") {
        if (mpi.rank==0) MPI_File_write_at(file, offset-delim.length(), end.c_str(), end.length(), MPI_CHAR, MPI_STATUS_IGNORE);
        offset += end.length() - delim.length();
    }
    MPI_File_seek_shared(file, offset, MPI_SEEK_CUR);
}

#else
template<class TLattice, typename T>
void writeArrayTxt(std::ofstream& file, std::vector<T> data, std::string fmt, std::string delim="\n", std::string end="NONE", int ndim=1, std::string dimDelim=" ") {
    // Convert data to plaintext
    std::string dataString = "";
    for (int iData=0; iData<TLattice::N; iData++) {
        dataString += formatPoint(&data[iData*ndim], fmt, dimDelim, delim, ndim);
    }
    file << dataString;
}
#endif


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
    sprintf(header, "# vtk DataFile Version 3.0\nSimulation data at timestep %d\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA %d\n", timestep, TLattice::LZ, TLattice::LY, TLattice::LX, nPoints);
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
        // Write data
        writeArrayTxt<TLattice>(file, data, "%13.8f", "\n", "NONE", directions, " ");
    } (), ...);
}


template<class TLattice>
template<typename... TParameter>
void SaveHandler<TLattice>::saveDAT(int timestep, TParameter&... params) {
    // File layout:
    // TITLE = [NAME]
    // VARIABLES = x, y, z, [param]...
    // ZONE  t= "solid", f= point, I= [LX], J= [LY], K= LZ
    // [x]\t[y]\t[z]\t[param]...
    // [x]\t[y]\t[z]\t[param]...
    // ...

    std::string filePrefix = "data";
    char filename[512];
    sprintf(filename, "%s/%s_%d.dat", mDataDir.c_str(), filePrefix.c_str(), timestep);
    #ifdef MPIPARALLEL
    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    #else
    std::ofstream file(filename);
    #endif

    // Parameter data
    std::string parameterHeader = "";
    int nPoints = TLattice::LXdiv * TLattice::LYdiv * TLattice::LZdiv;
    int nData = 3 * nPoints;
    ([&] {
        int directions = params.mDirections;
        nData += nPoints * directions;
        if (directions == 1) {
            parameterHeader += std::string(", ") + TParameter::mName;
        } else {
            for (int idir=0; idir<directions; idir++) {
                parameterHeader += ", " + (TParameter::mName + std::to_string(idir));
            }
        }
    } (), ...);

    // Write the header
    char text[1024];
    sprintf(text, "TITLE = Simulation data at timestep %d\nVARIABLES = x, y, z%s", timestep, parameterHeader.c_str());
    writeText(file, text);
    sprintf(text, "\nZONE  t=\"solid\", f= point, I= %d, J= %d, K= %d\n", TLattice::LX, TLattice::LY, TLattice::LZ);
    writeText(file, text);

    // Write the data
    std::vector<std::string> data;
    data.reserve(nData);
    int nWidth = std::to_string(std::max({TLattice::LX, TLattice::LY, TLattice::LZ})).length();
    for (int x=0; x<TLattice::subArray[0]; x++) {
        for (int y=0; y<TLattice::subArray[1]; y++) {
            for (int z=0; z<TLattice::subArray[2]; z++) {
                int k = computeK<TLattice>(x, y, z);
                std::stringstream pointData;
                // Global coordinates
                pointData << std::setw(nWidth) << x+TLattice::LXMPIOffset-TLattice::HaloXWidth << "\t";
                pointData << std::setw(nWidth) << y+TLattice::LYMPIOffset-TLattice::HaloYWidth << "\t";
                pointData << std::setw(nWidth) << z+TLattice::LZMPIOffset-TLattice::HaloZWidth;
                // Parameters
                ([&] {
                    int directions = params.mDirections;
                    for (int idir=0; idir<directions; idir++) {
                        auto value = params.getParameter()[k*directions+idir];
                        if (mMaskSolid && isMasked<TLattice>(k)) value = mMaskValue; // Apply the solid mask, if using
                        pointData << "\t" << std::setw(12) << value;
                    }
                } (), ...);
                data.push_back(pointData.str());
            }
        }
    }
    writeArrayTxt<TLattice>(file, data, "%s");
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
template<class TParameter, int TNumDir>
void SaveHandler<TLattice>::saveParameter(int timestep) {
    // Setup array for saving
    std::vector<typename TParameter::ParamType> param = TParameter::template get<TLattice,TNumDir>();
    if (mMaskSolid) applyMask<TLattice>(param, TNumDir, mMaskValue);

    // File
    char fdump[512];
    sprintf(fdump, "%s/%s_t%i.mat", mDataDir.c_str(), TParameter::mName, timestep);
    int fileOffset = sizeof(typename TParameter::ParamType) * mpi.rank * TNumDir * (TLattice::LX * TLattice::LY * TLattice::LZ) / mpi.size;

    #ifdef MPIPARALLEL
    // Parallel saving with MPI
    MPI_File fh;
    MPI_File_open(MPI_COMM_SELF, fdump, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_seek(fh, fileOffset, MPI_SEEK_SET);
    MPI_File_write(fh,&param[TLattice::HaloSize*TNumDir*TParameter::instances], TNumDir*(TLattice::N-2*TLattice::HaloSize)*TParameter::instances, MPI_DOUBLE, MPI_STATUSES_IGNORE);
    MPI_File_close(&fh);

    #else
    // Saving without MPI
    std::ofstream fs(fdump, std::ios::out | std::ios::binary);
    fs.seekp(fileOffset);
    for (int k = TLattice::HaloSize; k < TLattice::N-TLattice::HaloSize; k++) {
        for(int idx = 0; idx < TParameter::instances; idx++) {
            for(int idx2 = 0; idx2 < TNumDir; idx2++) {
                fs.write((char *)(&param[k*TNumDir*TParameter::instances+idx*TNumDir+idx2]), sizeof(typename TParameter::ParamType));
            }
        }
    };
    fs.close();
    #endif

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
