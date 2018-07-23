#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <time.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "mpi.h"
#include <omp.h>


#define MASTER						0
#define MAX_OMP					8
#define CUDA_BLOCK_SIZE		1024
#define CUDA_PART					10000
#define INPUT_FILE_NAME		"data.txt"
#define OUTPUT_FILE_NAME		"outPut.txt"


struct point_t
{
	double	x, y; 
	double	vX, vY; 
	int		cluster; 

}typedef point_t;

struct Cluster
{
	double	x, y;
	double	sumX, sumY;
	int		numOfPoints;

}typedef cluster_t;


// ------------------ cuda ------------------
cudaError_t movePointsCuda(point_t* points, int partSize, double dT);
cudaError_t findAllDistancesWithCuda(int size, point_t *points, double *radius);

// ------------------------ mpi data types ------------------------
MPI_Datatype createPointDataType();
MPI_Datatype createClusterDataType();

// ------------------------ send & recv with mpi ------------------------
void sendInitialDataToSlaves(int K, double LIMIT, int N, double T, double QM, double dT, int numprocs, int partSize, point_t *points, cluster_t *clusters, MPI_Datatype MPI_POINT, MPI_Datatype MPI_CLUSTER);
void recvInitialDataFromMaster(int *K, double *LIMIT, int *N, double *T, double *QM, double *dT, int *partSize, point_t **points, cluster_t **clusters, MPI_Datatype MPI_POINT, MPI_Datatype MPI_CLUSTER);
void recvUpdatesFromSlaves(int K, int numprocs, int partSize, point_t *points, cluster_t *clusters, int *flag_clusterChanged, MPI_Datatype MPI_POINT, MPI_Datatype MPI_CLUSTER);
void sendUpdatstoMaster(int K, int partSize, point_t *points, cluster_t *clusters, int *flag_clusterChanged, MPI_Datatype MPI_POINT, MPI_Datatype MPI_CLUSTER);

// ------------------------ memory ------------------------
void freeMemory(point_t * points, cluster_t *clusters);

// ------------------------ files ------------------------
void readFromFile(point_t **points, int *N, int *K, double *T, double *dT, double *LIMIT, double *QM);
void writeToFile(cluster_t *clusters, int K, double q, double t);

// ------------------------ kmeands ------------------------
void kmeans(double *t, int myid, int LIMIT, int T, int N, double * q, double QM, double dT, int K, int numprocs, int partSize, point_t *points, cluster_t *clusters, MPI_Datatype MPI_POINT, MPI_Datatype MPI_CLUSTER);
void initClusters(cluster_t **clusters, point_t *points, int K);
int getNearestClusterToPoint(point_t point, cluster_t *clusters, int K);
void clusterPoints(int K, int partSize, point_t *points, cluster_t *clusters, int *clusterChanged);
void updateCentroids(int myid, cluster_t *clusters, int K, int numprocs, MPI_Datatype MPI_CLUSTER);
void clearData(int K, cluster_t *clusters);
void getQ(int myid, double *q, point_t *points, int N, int K, int numprocs, cluster_t *clusters);
void findDiameters(point_t *points, int N, double *diameters, int K);

double getDistance(double x1, double x2, double y1, double y2);


// ------------------------ backup ------------------------
//void movePoints(int partSize, double dT, point_t *points);

//void findAllDistances(point_t *points, int N, double *diameters, int K, double *distances);

//void movePoints(int partSize, double dT, point_t *points)
//{
//	int i;
//	double newX, newY;
//
//	omp_set_num_threads(60); // 60 is the max number
//
//#pragma omp parallel for num_threads(60)
//	for (i = 0; i < partSize; i++)
//	{
//		points[i].x += points[i].vX*dT;
//		points[i].y += points[i].vY*dT;
//	}
//}

//void findAllDistances(point_t *points, int N, double *diameters, int K, double *distances)
//{
//	double		dist = 0;
//	int			cluster, i, j;
//
//	omp_set_num_threads(MAX_OMP);
//
//	// find for each point the distance of the farest point in the same cluster ---> O ( N * N / MAX_OMP) ?logn?
//
//#pragma omp parallel num_threads(MAX_OMP) private (i, j, cluster, dist)
//	{
//#pragma omp for 	
//		for (i = 0; i < N - 1; i++)
//		{
//			cluster = points[i].cluster;
//		
//			for (j = i + 1 ; j < N; j++)
//			{
//				if (points[j].cluster == cluster) //if points of the same cluster
//				{
//					getDistance(&dist, points[i].x, points[j].x, points[i].y, points[j].y);
//					if (dist > distances[i])
//						distances[i] = dist;
//				}
//			}
//		}
//	}
//}
//void getDiameters(double *diameters, point_t *points, int N, int K)
//{
//	double		dist = 0;
//	int			cluster, i, j;
//
//
//	for (i = 0; i <N - 1; i++)
//	{
//		cluster = points[i].cluster;
//
//		for (j = i + 1; j < N; j++)
//		{
//			if (points[j].cluster == cluster) //if points of the same cluster
//			{
//				getDistance(&dist, points[i].x, points[j].x, points[i].y, points[j].y);
//				if (dist > diameters[cluster])
//					diameters[cluster] = dist;
//			}
//		}
//	}
//}

