#include "Header.h"

int main(int argc, char* argv[])
{
	double start, QM, dT,T, q = DBL_MAX, LIMIT,t;
	int N, K, remainder, n = 0;
	point_t *points = NULL;
	cluster_t *clusters = NULL;
	
	//mpi 
	int numprocs, myid, partSize;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Status status;
	MPI_Datatype MPI_POINT;
	MPI_Datatype MPI_CLUSTER;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_POINT = createPointDataType();
	MPI_CLUSTER = createClusterDataType();
	start = MPI_Wtime();

	// ------------------ init all data ------------------
	if (myid == MASTER)
	{
		readFromFile(&points, &N, &K, &T, &dT, &LIMIT, &QM);
		initClusters(&clusters, points, K);
		partSize = N / numprocs;
		remainder = N % numprocs;
		
		sendInitialDataToSlaves(K, LIMIT, N, T, QM, dT, numprocs, partSize, points, clusters, MPI_POINT, MPI_CLUSTER);
	}
	else
		recvInitialDataFromMaster(&K, &LIMIT, &N, &T, &QM, &dT, &partSize, &points, &clusters, MPI_POINT, MPI_CLUSTER);


	// ------------------ start iteration in order to find q < QM ------------------
	do {
	
		//if its not first iteration we need to move the points
		if (myid == MASTER && n > 0)
			movePointsCuda(points + (numprocs - 1)*partSize, partSize + remainder , dT);  //  O( 1 ) because of cuda

		else if (n > 0) //slave 
			movePointsCuda(points, partSize, dT);	 //  O( 1 ) because of cuda									

		//movePoints(partSize, dT, points);							

		kmeans(&t,myid, LIMIT, T, N, &q, QM, dT, K, numprocs, partSize, points, clusters, MPI_POINT, MPI_CLUSTER); 

		getQ(myid, &q, points, N, K, numprocs, clusters);  //  O( N / CUDA_PART ) * O(K*K) 
		n++;
	} while ( (n*dT) < T && q > QM && n<3 );


	// ------------------ write output and finalize ------------------
	if (myid == MASTER)
	{
		writeToFile(clusters, K, q, (n*dT));
		printf("\n q: %lf\n t: %lf\n total time : %lf\n", q, n*dT, MPI_Wtime() - start);
	}
	freeMemory(points, clusters);
	MPI_Finalize();
	return 0;
}

void kmeans(double *t,int myid, int LIMIT, int T, int N, double * q, double QM, double dT, int K, int numprocs, int partSize, point_t *points, cluster_t *clusters, MPI_Datatype MPI_POINT, MPI_Datatype MPI_CLUSTER)
{
	int itr = 0;
	int remainder = N % numprocs;
	int flag_clusterChanged;	//flag == 1 if points moved between the clusters
	cluster_t *clustersFromSlave = (cluster_t*)calloc(K, sizeof(cluster_t));
	MPI_Status status;
	MPI_Comm comm = MPI_COMM_WORLD;
	cudaError_t cudaStatus;

	// cluster points and update centroids to be the avarege
	// iterate while points move between clusters && itr<LIMIT

	do {
		flag_clusterChanged = 0;
		clearData(K, clusters);

		if (myid == MASTER)
		{
			clusterPoints(K, partSize + remainder, points + (numprocs - 1)*partSize, clusters, &flag_clusterChanged);	//  O( partSize / MAX_OMP )
			recvUpdatesFromSlaves(K, numprocs, partSize, points, clusters, &flag_clusterChanged, MPI_POINT, MPI_CLUSTER);
		}
		else
		{
			clusterPoints(K, partSize, points, clusters, &flag_clusterChanged);											// ---> O( partSize / MAX_OMP )
			sendUpdatstoMaster(K, partSize, points, clusters, &flag_clusterChanged, MPI_POINT, MPI_CLUSTER);
		}

		updateCentroids(myid, clusters, K, numprocs, MPI_CLUSTER);														// ---> O ( K )
		itr++;
	} while (flag_clusterChanged > 0 && itr < LIMIT);
	(*t) = itr * dT;
}

double getDistance(double x1, double x2, double y1, double y2)
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}


// ------------------------ find quality ------------------------

void getQ(int myid, double *q, point_t *points, int N, int K, int numprocs, cluster_t *clusters)
{
	int			slave, i, j;
	double		*diameters, dist;
	MPI_Comm	comm = MPI_COMM_WORLD;
	MPI_Status	status;

	if (myid == MASTER)
	{
		// master calcs q and sends it to slaves
		// O( N / CUDA_PART ) * O(K*K)

		*q = 0;
		diameters = (double *)calloc(K, sizeof(double));
		findDiameters(points, N, diameters, K);	//  O( N / CUDA_PART )
		
		// O(k*k)
		for (i = 0; i < K; i++)
		{
			for (j = 0; j < K; j++)
			{
				if (i != j)
				{
					dist = getDistance(clusters[i].x, clusters[j].x, clusters[i].y, clusters[j].y);
					*q += diameters[i] / dist;
				}
			}
		}
		*q = *q / (K*(K - 1));

		for (slave = 1; slave < numprocs; slave++)
			MPI_Send(q, 1, MPI_DOUBLE, slave, 0, MPI_COMM_WORLD); //send q to slaves
	}
	else //slaves
		MPI_Recv(q, 1, MPI_DOUBLE, MASTER, 0, comm, &status);	 //resv q from master
}

void findDiameters(point_t *points, int N, double *diameters, int K)
{
	int			i;
	int			count = N / CUDA_PART;
	int			remainder = N % CUDA_PART;
	double		*distances = (double*)calloc(N, sizeof(double));

	//find max distances from each point in its class
	// send each time CUDA_PART so cudas memory not over-loaded

	for (i = 0; i < count; i++)		//	O( 1 )*O( count )
		findAllDistancesWithCuda( CUDA_PART, points + i*CUDA_PART, distances + i*CUDA_PART);
	if( remainder > 0 )
		findAllDistancesWithCuda( remainder, points + i*CUDA_PART, distances + i*CUDA_PART);

	// find max distance for each cluster
	for (i = 0; i < N; i++)	// O ( N )
	{
		if (distances[i] > diameters[points[i].cluster])
			diameters[points[i].cluster] = distances[i];
	}
}


// ------------------------ clusters related ------------------------

void clusterPoints(int K, int partSize, point_t *points, cluster_t *clusters, int *clusterChanged)
{
	//clustering points with omp 

	int			i, j, cluster, tid;
	int			*tempFlag = (int*)calloc(MAX_OMP, sizeof(int));
	cluster_t	*tempClusters = (cluster_t*)calloc(K * MAX_OMP, sizeof(cluster_t));

	omp_set_num_threads(MAX_OMP);

	// calc data into tempClusters ---> O( partSize / MAX_OMP )
#pragma omp parallel num_threads(MAX_OMP) private (cluster, tid)
	{
		tid = omp_get_thread_num();
#pragma omp for
		for (i = 0; i < partSize; i++)
		{
			cluster = getNearestClusterToPoint(points[i], clusters, K); 

			if (points[i].cluster != cluster)
				tempFlag[tid] = 1;	// points claster changed 	

			points[i].cluster = cluster;
			tempClusters[cluster + K*tid].numOfPoints++;
			tempClusters[cluster + K*tid].sumX += points[i].x;
			tempClusters[cluster + K*tid].sumY += points[i].y;
		}
	}

	// combine all data from tempClusters to clusters ---> O(K*MAX_OMP)
	for (i = 0; i < K; i++)
	{
		for (j = 0; j < MAX_OMP; j++)
		{
			clusters[i].numOfPoints += tempClusters[i + j*K].numOfPoints;
			clusters[i].sumX += tempClusters[i + j*K].sumX;
			clusters[i].sumY += tempClusters[i + j*K].sumY;
		}
	}
	for (i = 0; i < MAX_OMP; i++)
		(*clusterChanged) += tempFlag[i];
}

void initClusters(cluster_t **clusters, point_t *points, int K)
{
	int i;
	*clusters = (cluster_t*)calloc(K, sizeof(cluster_t));

	//init centroids with first k points and other data with 0 - O(K)
	for (i = 0; i < K; i++)
	{
		(*clusters)[i].x = points[i].x; 
		(*clusters)[i].y = points[i].y;
		(*clusters)[i].sumX = 0;
		(*clusters)[i].sumY = 0;
		(*clusters)[i].numOfPoints = 0;
	}
}

int getNearestClusterToPoint(point_t point, cluster_t *clusters, int K)
{
	int		i, nearestCluster = 0;
	double	dist, minDistance = DBL_MAX;

#pragma omp parallel for // O(1)
	for (i = 0; i < K; i++)
	{
		dist = getDistance(point.x, clusters[i].x, point.y, clusters[i].y );
		if (dist < minDistance)
		{
			nearestCluster = i;
			minDistance = dist;
		}
	}
	return nearestCluster;
}

void updateCentroids(int myid, cluster_t *clusters, int K, int numprocs, MPI_Datatype MPI_CLUSTER)
{
	int			slave, i, size;
	MPI_Comm	comm = MPI_COMM_WORLD;
	MPI_Status	status;

	//only master has all the clusters, he updates them and sends to slaves  - O( K )
	if (myid == MASTER)
	{
		for (i = 0; i < K; i++)
		{
			size = clusters[i].numOfPoints;
			if (size > 0)
			{
				clusters[i].x = clusters[i].sumX / size;
				clusters[i].y = clusters[i].sumY / size;
			}
		}

		//send new clusters to slaves 
		for (slave = 1; slave < numprocs; slave++)
			MPI_Send(clusters, K, MPI_CLUSTER, slave, 0, MPI_COMM_WORLD);
	}
	else
		MPI_Recv(clusters, K, MPI_CLUSTER, MASTER, 0, comm, &status);
}

void clearData(int K, cluster_t *clusters)
{
	int i;

	// O(K)
	for (i = 0; i < K; i++)
	{
		clusters[i].sumX = 0;
		clusters[i].sumY = 0;
		clusters[i].numOfPoints = 0;
	}
}


// ------------------------ send & recv with mpi ------------------------

void sendInitialDataToSlaves(int K, double LIMIT, int N, double T, double QM, double dT, int numprocs, int partSize, point_t *points, cluster_t *clusters, MPI_Datatype MPI_POINT, MPI_Datatype MPI_CLUSTER)
{
	int slave;

	for (slave = 1; slave < numprocs; slave++)
	{
		MPI_Send(&T, 1, MPI_DOUBLE, slave, 0, MPI_COMM_WORLD);
		MPI_Send(&LIMIT, 1, MPI_DOUBLE, slave, 0, MPI_COMM_WORLD);
		MPI_Send(&N, 1, MPI_INT, slave, 0, MPI_COMM_WORLD);
		MPI_Send(&QM, 1, MPI_DOUBLE, slave, 0, MPI_COMM_WORLD);
		MPI_Send(&dT, 1, MPI_DOUBLE, slave, 0, MPI_COMM_WORLD);
		MPI_Send(&K, 1, MPI_INT, slave, 0, MPI_COMM_WORLD);
		MPI_Send(&partSize, 1, MPI_INT, slave, 0, MPI_COMM_WORLD);

		MPI_Send(points + partSize*(slave - 1), partSize, MPI_POINT, slave, 0, MPI_COMM_WORLD);
		MPI_Send(clusters, K, MPI_CLUSTER, slave, 0, MPI_COMM_WORLD);
	}
}

void recvInitialDataFromMaster(int *K, double *LIMIT, int *N, double *T, double *QM, double *dT, int *partSize, point_t **points, cluster_t **clusters, MPI_Datatype MPI_POINT, MPI_Datatype MPI_CLUSTER)
{
	MPI_Comm		comm = MPI_COMM_WORLD;
	MPI_Status		status;

	//recv sizes
	MPI_Recv(T, 1, MPI_DOUBLE, MASTER, 0, comm, &status);
	MPI_Recv(LIMIT, 1, MPI_DOUBLE, MASTER, 0, comm, &status);
	MPI_Recv(N, 1, MPI_INT, MASTER, 0, comm, &status);
	MPI_Recv(QM, 1, MPI_DOUBLE, MASTER, 0, comm, &status);
	MPI_Recv(dT, 1, MPI_DOUBLE, MASTER, 0, comm, &status);
	MPI_Recv(K, 1, MPI_INT, MASTER, 0, comm, &status);
	MPI_Recv(partSize, 1, MPI_INT, MASTER, 0, comm, &status);

	//allocate memory
	*points = (point_t*)calloc(*partSize, sizeof(point_t));
	*clusters = (cluster_t*)calloc(*K, sizeof(cluster_t));

	//recv points and clusters
	MPI_Recv(*points, *partSize, MPI_POINT, MASTER, 0, comm, &status);
	MPI_Recv(*clusters, *K, MPI_CLUSTER, MASTER, 0, comm, &status);
}

void recvUpdatesFromSlaves(int K, int numprocs, int partSize, point_t *points, cluster_t *clusters, int *flag_clusterChanged, MPI_Datatype MPI_POINT, MPI_Datatype MPI_CLUSTER)
{
	int			i, j, flag_FromSlave;
	cluster_t	*clustersFromSlave = (cluster_t*)calloc(K, sizeof(cluster_t));
	MPI_Comm	comm = MPI_COMM_WORLD;
	MPI_Status	status;

	//master recvs data from slaves and updates the clusters 
	for (i = 1; i < numprocs; i++)
	{
		MPI_Recv(points + partSize*(i - 1), partSize, MPI_POINT, i, 0, comm, &status);//////any source?
		MPI_Recv(clustersFromSlave, K, MPI_CLUSTER, i, 0, comm, &status);
		MPI_Recv(&flag_FromSlave, 1, MPI_INT, i, 0, comm, &status);
		(*flag_clusterChanged) += flag_FromSlave;

		for (j = 0; j < K; j++)
		{
			clusters[j].numOfPoints += clustersFromSlave[j].numOfPoints;
			clusters[j].sumX += clustersFromSlave[j].sumX;
			clusters[j].sumY += clustersFromSlave[j].sumY;
		}
	}
	//master sends updated clusters and flag to slaves
	for (i = 1; i < numprocs; i++)
	{
		MPI_Send(clusters, K, MPI_CLUSTER, i, 0, MPI_COMM_WORLD);
		MPI_Send(flag_clusterChanged, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
}

void sendUpdatstoMaster(int K, int partSize, point_t *points, cluster_t *clusters, int *flag_clusterChanged, MPI_Datatype MPI_POINT, MPI_Datatype MPI_CLUSTER)
{
	MPI_Status	status;
	MPI_Comm	comm = MPI_COMM_WORLD;

	//slave sends his part to master
	MPI_Send(points, partSize, MPI_POINT, MASTER, 0, MPI_COMM_WORLD);
	MPI_Send(clusters, K, MPI_CLUSTER, MASTER, 0, MPI_COMM_WORLD);
	MPI_Send(flag_clusterChanged, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);

	//slave recvs updated clusters and flag
	MPI_Recv(clusters, K, MPI_CLUSTER, MASTER, 0, comm, &status);
	MPI_Recv(flag_clusterChanged, 1, MPI_INT, MASTER, 0, comm, &status);
}


// ------------------------ mpi data types ------------------------

MPI_Datatype createPointDataType()
{
	point_t point;
	MPI_Datatype MPI_POINT;
	int blocklen[5] = { 1,1,1,1,1 };
	MPI_Datatype type[5] = { MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE,MPI_INT };
	MPI_Aint disp[5] = { (char *)&point.x - (char *)&point , (char *)&point.y - (char *)&point,
						(char *)&point.vX - (char *)&point, (char *)&point.vY - (char *)&point,
						(char *)&point.cluster - (char *)&point };
	MPI_Type_create_struct(5, blocklen, disp, type, &MPI_POINT);
	MPI_Type_commit(&MPI_POINT);
	return MPI_POINT;
}

MPI_Datatype createClusterDataType()
{
	cluster_t cluster;
	MPI_Datatype MPI_CLUSTER;
	int blocklen[5] = { 1,1,1,1,1 };
	MPI_Datatype type[5] = { MPI_DOUBLE,MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE,MPI_INT };
	MPI_Aint disp[5] = { (char *)&cluster.x - (char *)&cluster , (char *)&cluster.y - (char *)&cluster,
						(char *)&cluster.sumX - (char *)&cluster, (char *)&cluster.sumY - (char *)&cluster,
						(char *)&cluster.numOfPoints - (char *)&cluster };
	MPI_Type_create_struct(5, blocklen, disp, type, &MPI_CLUSTER);
	MPI_Type_commit(&MPI_CLUSTER);
	return MPI_CLUSTER;
}


// ------------------------ files ------------------------

void readFromFile(point_t **points, int *N, int *K, double *T, double *dT, double *LIMIT, double *QM)
{
	int i;
	FILE* f = fopen(INPUT_FILE_NAME, "r");

	fscanf(f, "%d %d %lf %lf %lf %lf", N, K, T, dT, LIMIT, QM);

	*points = (point_t*)calloc(*N, sizeof(point_t));

	for (i = 0; i < *N; i++)
	{
		fscanf_s(f, "%lf %lf %lf %lf", &((*points)[i].x), &((*points)[i].y), &((*points)[i].vX), &((*points)[i].vY));
		(*points)[i].cluster = -1;
	}
	fclose(f);
}

void writeToFile(cluster_t *clusters, int K, double q, double t)
{
	int i;
	FILE* f = fopen(OUTPUT_FILE_NAME, "w");
	fprintf(f, "First occurrence at t = %lf with q =%lf \nCenters of the clusters:\n", t, q);

	for (i = 0; i < K; i++)
		fprintf(f, "%lf %f \n", clusters[i].x, clusters[i].y);

	fclose(f);
}


// ------------------------ memory ------------------------

void freeMemory(point_t * points, cluster_t *clusters)
{
	free(points);
	free(clusters);
}

