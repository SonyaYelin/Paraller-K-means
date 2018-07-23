#pragma once
#include "Header.h"

void free(point_t *point);
void free(point_t *point, double *distances);

__global__ void movePointsKernel(point_t* points, int partSize, double dT)
{
	int tid = threadIdx.x;
	int bid = blockIdx.x;
	int i = bid*CUDA_BLOCK_SIZE + tid;
	if (i < partSize)
	{
		points[i].x += points[i].vX*dT;
		points[i].y += points[i].vY*dT;
	}
}

__global__ void findAllDistances(int partSize, point_t *points, double *dDistances)
{
	int i;
	int tid = threadIdx.x;
	int bid = blockIdx.x;
	int index = bid*CUDA_BLOCK_SIZE + tid;

	double dist = 0, max = 0;

	//each thread gets a point and finds the distance to the point that is farest and has same cluster

	for (i = index + 1; i < partSize; i++)
	{
		if (points[index].cluster == points[i].cluster)
		{
			dist = sqrt(pow(points[index].x - points[i].x, 2) + pow(points[index].y - points[i].y, 2));
			if (dist > max)
				max = dist;
		}
	}

	if (index < partSize)
		dDistances[index] = max;
}

cudaError_t movePointsCuda(point_t* points, int partSize, double dT)
{
	int numOfBlocks;
	point_t *dPoints;
	cudaError_t cudaStatus;
	numOfBlocks = partSize / CUDA_BLOCK_SIZE;
	if (partSize%CUDA_BLOCK_SIZE > 0)
		numOfBlocks += 1;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		free(dPoints);
	}

	// Allocate GPU buffers   .
	cudaStatus = cudaMalloc((void**)&dPoints, partSize * sizeof(point_t));
	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, "cudaMalloc dPoints failed!");
		free(dPoints);
	}
	// Copy input vectors from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dPoints, points, partSize * sizeof(point_t), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, "cudaMemcpy dPoints failed!");
		free(dPoints);
	}

	// Launch a kernel on the GPU with one thread for each element.
	movePointsKernel << <numOfBlocks, CUDA_BLOCK_SIZE >> >(dPoints, partSize, dT);
	
	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "movePoints launch failed: %s\n", cudaGetErrorString(cudaStatus));
		free(dPoints);
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching movePoints!\n", cudaStatus);
		free(dPoints);
	}

	// Copy output vector from GPU buffer to host memory.
	//printf("\n(%f,%f)\n", dev_points[0].x, dev_points[0].y);
	cudaStatus = cudaMemcpy(points, dPoints, partSize * sizeof(point_t), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, "cudaMemcpy from device  failed!");
		free(dPoints);
	}

	free(dPoints);
	return cudaStatus;
}

cudaError_t findAllDistancesWithCuda(int partSize, point_t *points, double *distances)
{
	cudaError_t	cudaStatus;
	point_t		*dPoints;
	double		*dDistances;
	int			numOfBlocks = partSize / CUDA_BLOCK_SIZE;
	if (partSize % CUDA_BLOCK_SIZE > 0)
		numOfBlocks++;

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		free(dPoints, dDistances);
	}

	// Allocate GPU buffers 
	cudaStatus = cudaMalloc((void**)&dPoints, partSize * sizeof(point_t));
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaMalloc dPoints failed!");
		free(dPoints, dDistances);
	}
	cudaStatus = cudaMalloc((void**)&dDistances, partSize * sizeof(double));
	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, "cudaMalloc dDistances failed!");
		free(dPoints, dDistances);
	}
	// Copy input from host memory to GPU buffers.
	cudaStatus = cudaMemcpy(dPoints, points, partSize * sizeof(point_t), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy dPoints from host failed!");
		free(dPoints, dDistances);
	}

	cudaStatus = cudaMemcpy(dDistances, distances, partSize * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, "cudaMemcpy dDistances  from host failed!");
		free(dPoints, dDistances);
	}

	// Launch a kernel on the GPU with one thread for each element.
	findAllDistances << <numOfBlocks, CUDA_BLOCK_SIZE >> >(partSize, dPoints, dDistances);

	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d!\n", cudaStatus);
		free(dPoints, dDistances);
	}

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, " launch failed: %s\n", cudaGetErrorString(cudaStatus));
		free(dPoints, dDistances);
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d!\n", cudaStatus);
		free(dPoints, dDistances);
	}

	// Copy output vector from GPU buffer to host memory.
	//printf("\n(%f,%f)\n", dev_points[0].x, dev_points[0].y);
	cudaStatus = cudaMemcpy(distances, dDistances, partSize * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) 
	{
		fprintf(stderr, "cudaMemcpy from device  failed!");
		free(dPoints, dDistances);
	}
		
	free(dPoints, dDistances);
	return cudaStatus;
}

void free(point_t * point)
{
	cudaFree(point);
}

void free(point_t * point, double *distances)
{
	cudaFree(point);
	cudaFree(distances);
}