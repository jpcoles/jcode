/*
 * Perform memory verification and transfer tests on graphics cards supporting 
 * the nVidia CUDA API.
 *
 * Written by Jonathan Coles
 *
 * Build with: make clean && make
 *
 * Last Modified On: 21-DEC-2007
 */

#include <stdio.h>
#include <cuda.h>
#include <assert.h>
#include <cutil.h>
#include "timing.h"

CPUDEFS

#define CHECK(err, msg) if (ret != cudaSuccess) { printf msg; goto error; }

cudaError_t MALLOC_HOST(void **ptr, size_t size)
{
    *ptr = malloc(size);
    return *ptr != NULL ? cudaSuccess : cudaErrorMemoryAllocation;
}

cudaError_t MEMSET_HOST(void *ptr, int val, size_t size) 
{ 
    memset(ptr, val, size); 
    return cudaSuccess;
}
cudaError_t FREE_HOST(void *ptr)
{ 
    free(ptr);
    return cudaSuccess;
}

int verify_copy_tests(struct cudaDeviceProp *prop)
{
    long i, j, nbytes;
    unsigned char *host = NULL, *dev = NULL;
    cudaError_t ret = cudaSuccess;
    
    const long ntransfers = 3;
    const long max_nbytes = 2*prop->totalGlobalMem - 1;

    for (nbytes=1; nbytes <= max_nbytes; nbytes *= 2)
    {
        if (nbytes > prop->totalGlobalMem) nbytes = prop->totalGlobalMem;

        printf("Verifying %i byte transfers...\n", nbytes);
        ret = MALLOC_HOST((void **)&host, nbytes); 
        CHECK(ret, ("Can't allocate %i bytes on host\n", nbytes));
        ret = MEMSET_HOST(host, (nbytes & 0xFF), nbytes);
        CHECK(ret, ("Can't set %i bytes on host\n", nbytes));
        ret = cudaMalloc((void **)&dev, nbytes);
        CHECK(ret, ("Can't allocate %i bytes on device.\n", nbytes));
        for (j=0; j < ntransfers; j++)
        {
            ret = cudaMemcpy(dev, host, nbytes, cudaMemcpyHostToDevice);
            CHECK(ret, ("Can't copy %i bytes to device\n", nbytes));
            MEMSET_HOST(host, 0, nbytes); /* In case the memcpy is just failing */
            ret = cudaMemcpy(host, dev, nbytes, cudaMemcpyDeviceToHost);
            CHECK(ret, ("Can't copy %i bytes to host\n", nbytes));
            for (i=0; i < nbytes; i++) 
            {
                if (host[i] != (nbytes & 0xFF))
                {
                    printf("Mismatched bytes after transfer (%i: %i, %i)\n", 
                           i, host[i], (nbytes & 0xFF));
                    goto error;
                }
            }
        }
        ret = cudaFree(dev); dev = NULL; CHECK(ret, (""));
        ret = FREE_HOST(host); dev = NULL; CHECK(ret, (""));
    }

error:

    if (host != NULL) FREE_HOST(host);
    if (dev  != NULL) cudaFree(dev);

    if (ret != cudaSuccess) printf("CUDA Error: %s\n", cudaGetErrorString(ret));

    return ret != cudaSuccess;
}

int transfer_rate_tests(struct cudaDeviceProp *prop)
{
    long j, nbytes;
    unsigned char *host = NULL, *dev = NULL;
    double start, end;
    cudaError_t ret;
    
    const long ntransfers = 10000;
    const long max_nbytes = 2*prop->totalGlobalMem - 1;

    for (nbytes=1; nbytes <= max_nbytes; nbytes *= 2)
    {
        if (nbytes > prop->totalGlobalMem) nbytes = prop->totalGlobalMem;

        printf("Timing %i byte transfers... ", nbytes);
        ret = MALLOC_HOST((void **)&host, nbytes); 
        CHECK(ret, ("Can't allocate %i bytes on host\n", nbytes));
        ret = cudaMalloc((void **)&dev, nbytes);
        CHECK(ret, ("Can't allocate %i bytes on device.\n", nbytes));
        
        start = CPUTIME;
        for (j=0; j < ntransfers; j++)
            ret = cudaMemcpy(dev, host, nbytes, cudaMemcpyHostToDevice);
        end = CPUTIME;
        printf("%e ", (end-start) / ntransfers);
        printf("%f ", (double)nbytes / ((end-start) / (double)ntransfers) / (1024*1024*1024.0));

        start = CPUTIME;
        for (j=0; j < ntransfers; j++)
            ret = cudaMemcpy(host, dev, nbytes, cudaMemcpyDeviceToHost);
        end = CPUTIME;
        printf("%e ", (end-start) / ntransfers);
        printf("%f ", (double)nbytes / ((end-start) / (double)ntransfers) / (1024*1024*1024.0));

        printf("\n");

        cudaFree(dev); dev = NULL;
        FREE_HOST(host); host = NULL;
    }

error:

    if (host != NULL) FREE_HOST(host);
    if (dev  != NULL) cudaFree(dev);

    return 0;
}

int main(int argc, char **argv)
{
    int i;
    int deviceCount;
    struct cudaDeviceProp prop;
    cudaError_t err;

    cudaGetDeviceCount(&deviceCount);

    printf("Found %i device(s).\n", deviceCount);

    for (i=0; i < deviceCount; i++)
    {
        printf("--------------------------------------------\n"
               "Device %i\n", i);

        err = cudaGetDeviceProperties(&prop, i);
        if (err != cudaSuccess)
        {
            printf("Error %i\n", err);
            continue;
        }

        printf("Name: '%s'\n"
               "Total Global Memory: %u\n"
               "Shared Memory Per Block: %u\n"
               "Registers Per Block: %i\n"
               "Warp Size: %i\n"
               "Memory Pitch: %u\n"
               "Maximum Threads Per Block: %i\n"
               "Maximum Size of Each Block Dimension: %i %i %i\n"
               "Maximum Size of Each Grid Dimension: %i %i %i\n"
               "Total Constant Memory: %i\n"
               "Revision: %i.%i\n"
               "Clockrate: %iHz\n",
               prop.name,
               prop.totalGlobalMem,
               prop.sharedMemPerBlock,
               prop.regsPerBlock,
               prop.warpSize,
               prop.memPitch,
               prop.maxThreadsPerBlock,
               prop.maxThreadsDim[0],
               prop.maxThreadsDim[1],
               prop.maxThreadsDim[2],
               prop.maxGridSize[0],
               prop.maxGridSize[1],
               prop.maxGridSize[2],
               prop.totalConstMem,
               prop.major,
               prop.minor,
               prop.clockRate);
        printf("Size of Property Structure: %u\n", sizeof(prop));

        cudaSetDevice(i);
        verify_copy_tests(&prop);
        transfer_rate_tests(&prop);
    }

    return 0;
}

