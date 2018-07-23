//#include <stdlib.h>
//#include <stdio.h>
//#include <float.h>
//#include <math.h>
//#include <time.h>
//
//#define MASTER 0
//
//struct Point
//{
//	double	x, y;
//	double	vX, vY;
//	double	currentPosX, currentPosY;
//	int		cluster;
//
//}typedef point_t;
//
//struct Cluster
//{
//	double x, y;
//	double sumX, sumY;
//	int numOfPoints;
//
//}typedef cluster_t;
//
//void readFromFile(point_t **points, int *N, int *K, int *T, double *dT, int *LIMIT, double *QM)
//{
//	int i;
//	FILE* f = fopen("input.txt", "r");
//	fscanf(f, "%d %d %d %lf %d %lf", N, K, T, dT, LIMIT, QM);
//
//	*points = (point_t*)calloc(*N, sizeof(point_t));
//
//	for (i = 0; i < *N; i++)
//	{
//		fscanf_s(f, "%lf %lf %lf %lf", &((*points)[i].x), &((*points)[i].y), &((*points)[i].vX), &((*points)[i].vY));
//		(*points)[i].cluster = -1;
//	}
//	fclose(f);
//}
//
//void writeToFile(cluster_t *clusters, int K, double q, double t)
//{
//	int i;
//	FILE* f = fopen("outPut.txt", "w");
//
//	fprintf(f, "First occurrence at t = %lf with q =%lf \nCenters of the clusters:\n", t, q);
//
//	for (i = 0; i < K; i++)
//		fprintf(f, "%lf %f \n", clusters[i].x, clusters[i].y);
//
//	fclose(f);
//}
//
//void initClusters(cluster_t **clusters, point_t *points, int K)
//{
//	int i;
//	*clusters = (cluster_t*)calloc(K, sizeof(cluster_t));
//
//	for (i = 0; i < K; i++)
//	{
//		(*clusters)[i].x = points[i].x;
//		(*clusters)[i].y = points[i].y;
//		(*clusters)[i].sumX = 0;
//		(*clusters)[i].sumY = 0;
//		(*clusters)[i].numOfPoints = 0;
//	}
//}
//
//void movePoints(point_t *points, int size, double t)
//{
//	int i;
//
//	for (i = 0; i < size; i++)
//	{
//		points[i].currentPosX = points[i].x + t*points[i].vX; // x(t) = x0+v0*t 
//		points[i].currentPosY = points[i].y + t*points[i].vY;
//	}
//}
//
//void getDistance(double *dist, double x1, double x2, double y1, double y2)
//{
//	*dist = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
//}
//
//int getNearestClusterToPoint(point_t point, cluster_t *clusters, int K)
//{
//	int		i, nearestCluster = 0;
//	double	dist, minDistance = DBL_MAX;
//
//	for (i = 0; i < K; i++)
//	{
//		getDistance(&dist, point.currentPosX, clusters[i].x, point.currentPosY, clusters[i].y);
//		if (dist < minDistance)
//		{
//			nearestCluster = i;
//			minDistance = dist;
//		}
//	}
//	//printf(" \n min dist: %lf", minDistance);
//	return nearestCluster;
//}
//
//void clusterPoints(int K, int partSize, point_t *points, cluster_t *clusters, int *clusterChanged)
//{
//	int i, cluster;
//
//	printf("--------------------------------------------------------\n");
//	for (i = 0; i < partSize; i++)
//	{
//		cluster = getNearestClusterToPoint(points[i], clusters, K);
//
//		if (points[i].cluster != cluster)
//			(*clusterChanged) = 1;
//
//		points[i].cluster = cluster;
//		clusters[cluster].numOfPoints++;
//		clusters[cluster].sumX += points[i].currentPosX;
//		clusters[cluster].sumY += points[i].currentPosY;
//
//		printf("cur pos: (%lf, %lf) \n", points[i].currentPosX, points[i].currentPosY);
//		printf("sum x: %lf, sum y : %lf \n", clusters[cluster].sumX, clusters[cluster].sumY);
//
//	}
//	system("pause");
//}
//
//void updateCentroids(cluster_t *clusters, int K)
//{
//	int i, size;
//
//	for (i = 0; i < K; i++)
//	{
//		size = clusters[i].numOfPoints;
//		if (size > 0)
//		{
//			clusters[i].x = clusters[i].sumX / size;
//			clusters[i].y = clusters[i].sumY / size;
//		}
//		clusters[i].sumX = 0;
//		clusters[i].sumY = 0;
//		clusters[i].numOfPoints = 0;
//	}
//
//}
//
//void getDiameters(double *diameters, point_t *points, int N, int K)
//{
//	double		dist = 0;
//	int			cluster, i, j;
//
//	for (i = 0; i < N - 1; i++)
//	{
//		cluster = points[i].cluster;
//
//		for (j = i + 1; j < N; j++)
//		{
//			if (points[j].cluster == cluster) //if points of the same cluster
//			{
//				getDistance(&dist, points[i].currentPosX, points[j].currentPosX, points[i].currentPosY, points[j].currentPosY);
//				if (dist > diameters[cluster])
//					diameters[cluster] = dist;
//			}
//		}
//	}
//}
////
//double getQ(point_t *points, int N, int K, cluster_t *clusters)
//{
//	int			i, j, counter = 0;
//	double		dist, q = 0;
//	double		*diameters = (double *)calloc(K, sizeof(double));
//
//	getDiameters(diameters, points, N, K);
//
//	for (i = 0; i < K; i++)
//	{
//		for (j = i + 1; j < K; j++)
//		{
//			counter += 2;
//			getDistance(&dist, clusters[i].x, clusters[j].x, clusters[i].y, clusters[j].y);
//			q += (diameters[i] + diameters[j]) / dist;
//		}
//	}
//	return q / counter;
//}
//
//void freeMemory(point_t * points, cluster_t *clusters)
//{
//	free(points);
//	free(clusters);
//}
//
//void kmeans(int LIMIT, int T, int N, double *t, double *q, double QM, double dT, int K, point_t *points, cluster_t *clusters)
//{
//	int	flag_clusterChanged, itr, n = 0, i;
//	int cluster, size;
//	int c = 0;
//	*t = 0;
//
//	do
//	{
//		movePoints(points, N, *t);
//		itr = 0;
//		do
//		{
//
//			flag_clusterChanged = 0;
//			for (i = 0; i < K; i++)
//			{
//				clusters[i].sumX = 0;
//				clusters[i].sumY = 0;
//				clusters[i].numOfPoints = 0;
//			}
//			clusterPoints(K, N, points, clusters, &flag_clusterChanged); //if there is a point that moved to another cluster flag=1		
//
//			for (i = 0; i < N; i++)
//			{
//
//				cluster = getNearestClusterToPoint(points[i], clusters, K);
//
//				if (points[i].cluster != cluster)
//					(flag_clusterChanged) = 1;
//
//				points[i].cluster = cluster;
//				clusters[cluster].numOfPoints++;
//				clusters[cluster].sumX += points[i].currentPosX;
//				clusters[cluster].sumY += points[i].currentPosY;
//			}
//
//			updateCentroids(clusters, K);
//
//			itr++;
//
//		} while (flag_clusterChanged > 0 && itr < LIMIT);
//
//
//		(*q) = getQ(points, N, K, clusters);
//		n++;
//		*t = n*dT;
//		LIMIT--;
//	} while (*q > QM && *t < T);
//
//	printf("  \nq: %lf , t : %lf \n", *q, *t);
//}
//
//int main(int argc, char* argv[])
//{
//	double				QM, dT, t, q = DBL_MAX;
//	int					N, K, LIMIT, T;
//	point_t				*points = NULL;
//	cluster_t			*clusters = NULL;
//	clock_t				start = clock(), end;
//
//	readFromFile(&points, &N, &K, &T, &dT, &LIMIT, &QM);
//	initClusters(&clusters, points, K);
//
//	kmeans(LIMIT, T, N, &t, &q, QM, dT, K, points, clusters);
//
//	writeToFile(clusters, K, q, t);
//	freeMemory(points, clusters);
//	return 0;
//}
