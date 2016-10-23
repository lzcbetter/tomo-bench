#include "tomo_io.hpp"
#include <mpi.h>
//#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <stdio.h>
using namespace std;
/*
***************************************************************************************************
* func   name: env_init
* description: benchmark environment initialization 
*
* parameters : 
*             number of processes and rank
* return: none
***************************************************************************************************
*/
void env_init(int np, int rank, unsigned int P, unsigned int S, unsigned int C, unsigned int write_flag)
{
    char filename[100];
    #if _DEBUG_
        sprintf( filename, "Rank-%d-stdout.txt", rank );
        freopen( filename, "w", stdout );
    #else
        if(rank == 0)
        {
            sprintf( filename, "Rank-%d-NP-%d-P%d-S%d-C%d-TYPE%d.txt", rank, np, P, S, C, write_flag );
            //freopen( filename, "w", stdout );
        }
    #endif
}

/*
***************************************************************************************************
* func   name: tomo_io_finalize
* description: finalize, free allocated memory
*
* parameters : 
*             none
* return: none
***************************************************************************************************
*/
void tomo_io_finalize()
{
    //delete pin_buf;
    delete pout_buf;
}
/*
***************************************************************************************************
* func   name: memory_allocation
* description: allocate memory for hosting input/output
*
* parameters : 
*             # of projection, sinogram, column per rank
* return: none
***************************************************************************************************
*/
void memory_allocation(unsigned int p, unsigned int s, unsigned int c)
{
    unsigned long in_size  = s * p * c;               // input buffer size, one process
    unsigned long out_size = s * c * c;               // output buffer size, one process
    //pin_buf  = new float[in_size];
    pout_buf = new float[out_size];
}
/*
***************************************************************************************************
* func   name: main
* description: the main access function, it has several arguments. must be given in execution, 
*
* parameters : 
*             none
* return: none
***************************************************************************************************
*/
int main(int argc, char *argv[])
{
    int rank, np;
    MPI_File fh;
    int errcode;
    MPI_Offset offset;
    MPI_Status status;
    long long count;
    unsigned int write_flag;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(argc != 5){
        cout << "no sufficient parameters given" << endl;
        cout << "arguments should be given as: (MPIO write way, 0 for independent and 1 for collective) P, S, C" << endl;
        exit(-1);
    }else{
        write_flag = atoi(argv[1]);
        P = atoi(argv[2]);
        S = atoi(argv[3]);
        C = atoi(argv[4]);
    }
    if(S % np != 0){
        cout << "number of sinograms should be perfectly divisible by the number of processes" << endl;
        exit(-1);
    }
    env_init(np, rank, P, S, C, write_flag);
    memory_allocation(P, S/np, C);        // allocate memory to emulate output data access
    long long stride = sizeof(float) * C * C;
    offset = (long long)rank * (S/np) * stride;
    count = (S/np) * C * C;
    struct timeval fwrite_s, fwrite_e;
    double fwrite_eps;

    gettimeofday(&fwrite_s, NULL);                // start to count the time takes on writing
    char out_filename[256];
    sprintf( out_filename, "/projects/SDAV/zliu/tomo_out-NP-%d-P%d-S%d-C%d-TYPE%d", np, P, S, C, write_flag );
    //char *out_filename = (char *)"tomo_out";
    errcode = MPI_File_open(MPI_COMM_WORLD, out_filename,  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if (errcode != MPI_SUCCESS) {
		cout << "file open failed" << endl;
        exit(-1);
    }
    if (write_flag == 0){
        MPI_File_write_at(fh, offset, pout_buf, count, MPI_FLOAT, &status);
    }else{
        MPI_File_write_at_all(fh, offset, pout_buf, count, MPI_FLOAT, &status);
    }
    MPI_File_close(&fh);

    MPI_Barrier( MPI_COMM_WORLD );                // should wait untill all processes finish writing
    gettimeofday(&fwrite_e, NULL);                // stop count the time takes on writing
    fwrite_eps = fwrite_e.tv_sec - fwrite_s.tv_sec + (fwrite_e.tv_usec - fwrite_s.tv_usec) / 1e6;
    tomo_io_finalize();
    MPI_Finalize();
    if(rank == 0){
        remove(out_filename);                     // delete file to avoid wasting storage
        if (write_flag == 0){
            cout << "file write independently " << endl;
        }else{
            cout << "file write collectively " << endl;
        }
        cout << "P= " << P << ", S= " << S << ", C= " << C << endl;
        cout << "there are " << np << " processes, each process writes: " << count*sizeof(float) << " bytes" << endl; 
        cout << "total file size is: " << (float)np*count*sizeof(float)/1024.0/1024.0/1024.0 << " GB" << endl;
        cout << "file write time: " << fwrite_eps << " seconds" << endl;
    }
}