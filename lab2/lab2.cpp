#include <float.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <mpi.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <pthread.h>
#include <iomanip>

using namespace std;

#define NUM_THREADS 4

MPI_Datatype pointType;

struct Point {
	float coord[2];
	int index;
};

bool operator<(const Point& point1, const Point& point2) {
	return point1.coord[0] < point2.coord[0];
}

float getNextFloat() {
	return static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 1000000.0));
}

Point createFictivePoint() {
	Point fictivePoint;
	fictivePoint.coord[0] = FLT_MAX;
	fictivePoint.coord[1] = FLT_MAX;
	fictivePoint.index = -1;
	return fictivePoint;
}

void generatePoints(Point *points, int n, int i, int elementsNumber) {
	for(int j = 0; j < n; j++) {
		int index = n * i + j;
		if (index >= elementsNumber) {
			points[j] = createFictivePoint();
		}
		else {
			points[j].coord[0] = getNextFloat();
    			points[j].coord[1] = getNextFloat();
    			points[j].index = index;
		}
	}
}

void printPoints(Point *points, int n, int rank) {
	for(int i = 0; i < n; i++) {
    //if (points[i].index != -1) {
   	 cout << "x[" << rank <<  "][" << points[i].index << "] = " << fixed << setprecision(4) << points[i].coord[0] << endl;
    //}
	}
}

bool checkSorted(Point *points, int n) {
	for (int i = 0; i < n - 1; i++)
    	if (points[i+1] < points[i])
        	return false;
	return true;
}

void merge(const vector<int> &arr1, const vector<int> &arr2, vector<int> &result) {
	if (arr1.size() == 0 || arr2.size() == 0) return;
	else if (arr1.size() == 1 && arr2.size() == 1) {
    	result.push_back(arr1[0]);
    	result.push_back(arr2[0]);
    	//outputFile << ' ' << arr1[0] << ' ' << arr2[0] << endl;
    	//cout << ' ' << arr1[0] << ' ' << arr2[0] << endl;
    	//numberOfComparators++;
    	//int maxBeats = max(++lines[arr1[0]], ++lines[arr2[0]]);
    	//++lines[arr1[0]]; // = maxBeats;
    	//++lines[arr2[0]]; // = maxBeats;
    	return;
	}

	vector<int> odd1, odd2, even1, even2;

	for (int i = 0; i < arr1.size(); i++) {
    	if (i % 2 == 0) even1.push_back(arr1[i]);
    	else odd1.push_back(arr1[i]);
	}

	for (int i = 0; i < arr2.size(); i++) {
    	if (i % 2 == 0) even2.push_back(arr2[i]);
    	else odd2.push_back(arr2[i]);
	}

	merge(odd1, odd2, result);
	merge(even1, even2, result);

	vector<int> closing(arr1.begin(), arr1.end());
	closing.insert(closing.end(), arr2.begin(), arr2.end());

	for (int i = 1; i + 1 < closing.size(); i += 2) {
    	result.push_back(closing[i]);
    	result.push_back(closing[i + 1]);
    	//outputFile << ' ' << closing[i] << ' ' << closing[i + 1] << endl;
    	//cout << ' ' << closing[i] << ' ' << closing[i + 1] << endl;
    	//numberOfComparators++;
    	//int maxBeats = max(++lines[closing[i]], ++lines[closing[i + 1]]);
    	//++lines[closing[i]];// = maxBeats;
    	//++lines[closing[i + 1]];// = maxBeats;
	}
}

void batcherSort(const vector<int> &arr, vector<int> &result) {
	if (arr.size() == 1) return;
	vector<int> arr1(arr.begin(), arr.begin() + arr.size() / 2);
	vector<int> arr2(arr.begin() + arr.size() / 2, arr.end());
	batcherSort(arr1, result);
	batcherSort(arr2, result);
	merge(arr1, arr2, result);
}

vector<int> createSortingNet(int offset, int numprocs) {
	vector<int> sortingNet;
	vector<int> processorRanks;
   	for (int i = 0; i < numprocs; i++)
   		processorRanks.push_back(i + offset);
	sortingNet.clear();
   	batcherSort(processorRanks, sortingNet);
	return sortingNet;
}

void heapify(Point *points, int n, int i) {
	int max = i;
	int left = 2 * i + 1;
	int right = 2 * i + 2;

	if (left < n && points[max] < points[left])
    	max = left;

	if (right < n && points[max] < points[right])
    	max = right;

	if (max != i) {
    		Point temp = points[i];
    		points[i] = points[max];
    		points[max] = temp;
    		heapify(points, n, max);
	}
}

void buildHeap(Point *points, int n) {
	for (int i = n / 2; i >= 0; i--)
    	heapify(points, n, i);
}

void heapSort(Point *points, int n) {
	buildHeap(points, n);
		for (int i = n - 1; i >= 0; i--) {
    		Point temp = points[i];
    		points[i] = points[0];
    		points[0] = temp;
    		heapify(points, i, 0);
	}
}

void createPointDatatype() {
	Point point;
	int count = 2;
	MPI_Datatype datatypes[2] = {MPI_INT, MPI_INT};
	int lengths[2] = {2,1};
	MPI_Aint offsets[2];
	offsets[0] = offsetof(Point, coord);
	offsets[1] = offsetof(Point, index);
	MPI_Type_create_struct(count, lengths, offsets, datatypes, &pointType);
	MPI_Type_commit(&pointType);
}

void mergeSortedParts(Point *points, int l, int n1, int n2) {
	Point *left = new Point[n1], *right = new Point[n2];
	for (int i = l, j = 0; i < l + n1; i++, j++)
    	left[j] = points[i];
	for (int i = l + n1, j = 0; i < l + n1 + n2; i++, j++)
    	right[j] = points[i];
	for (int i = 0, j = 0, k = l; k < l + n1 + n2; k++) {
    	if (i >= n1) points[k] = right[j++];
    	else if (j >= n2) points[k] = left[i++];
    	else if (left[i] < right[j]) points[k] = left[i++];
    	else points[k] = right[j++];
	}
	delete [] left;
	delete [] right;
}

void mergeBlocks(Point *points, int numberOfBlocks, int *blockSizes) {
	if (numberOfBlocks == 1) return;
	int newNumberOfBlocks = ceil(numberOfBlocks / (float) 2);
	int newNumberOfNotFullBlocks = 0;
	int *newBlockSizes = new int[newNumberOfBlocks];
	int currentSize = 0;
	for (int i = 0, j = 0; i < numberOfBlocks; i += 2, j++) {
    	if (i + 1 < numberOfBlocks) {
        	newBlockSizes[j] = blockSizes[i] + blockSizes[i+1];
        	mergeSortedParts(
            	points, currentSize,
            	blockSizes[i], blockSizes[i+1]
        	);
        	currentSize = currentSize + blockSizes[i] + blockSizes[i+1];
    	}
    	else {
        	newBlockSizes[j] = blockSizes[i];
        	currentSize += blockSizes[i];
    	}
	}
	mergeBlocks(points, newNumberOfBlocks, newBlockSizes);
}

struct PthreadArgs {
	Point *points;
	int size;
};

void *sortArray(void *args) {
	PthreadArgs *pointsAndSize = (struct PthreadArgs *)args;
	sort(pointsAndSize->points, pointsAndSize->points + pointsAndSize->size);
	//heapSort(pointsAndSize->points, pointsAndSize->size);
	pthread_exit(NULL);
	return NULL;
}

void sortArrayBlocks(Point *points, int numberOfBlocks, int numberOfElementsPerFullBlock, int numberOfNotFullBlocks, int *blockSizes) {
	pthread_t *threads = new pthread_t[numberOfBlocks];
	int currentSize = 0;
	for (int i = 0; i < numberOfBlocks; i++) {
    	PthreadArgs *arrWithSize = new PthreadArgs();
    	if (i < numberOfBlocks - numberOfNotFullBlocks) {
        	arrWithSize->points = &points[currentSize];
        	currentSize += numberOfElementsPerFullBlock;
    	}
    	else {
        	arrWithSize->points = &points[currentSize];
        	currentSize += (numberOfElementsPerFullBlock - 1);
    	}
    	arrWithSize->size = blockSizes[i];
    	pthread_create(&threads[i], NULL, sortArray, (void *) arrWithSize);
	}
	for (int i = 0; i < numberOfBlocks; i++)
    	pthread_join(threads[i], NULL);
}

void parallelSort(Point *points, int n) {
	int numberOfBlocks = NUM_THREADS;
	int numberOfElementsPerFullBlock = ceil(n / (double) numberOfBlocks);
	int numberOfNotFullBlocks = numberOfBlocks * numberOfElementsPerFullBlock - n;

	int *blockSizes = new int[numberOfBlocks];
	for (int i = 0; i < numberOfBlocks; i++) {
    	if (i < numberOfBlocks - numberOfNotFullBlocks)
        	blockSizes[i] = numberOfElementsPerFullBlock;
    	else blockSizes[i] = numberOfElementsPerFullBlock - 1;
	}

	sortArrayBlocks(points, numberOfBlocks, numberOfElementsPerFullBlock, numberOfNotFullBlocks, blockSizes);
	mergeBlocks(points, numberOfBlocks, blockSizes);
}

void distributedSort(Point **points, int elementsNumber, int rank, vector<int> &sortingNet) {
	parallelSort(*points, elementsNumber);
	for (int i = 0; i < sortingNet.size(); i += 2) {
    	if (rank == sortingNet[i]) {
        	Point *newPoints = new Point[elementsNumber], *finalPoints = new Point[elementsNumber];
        	MPI_Send(*points, elementsNumber, pointType, sortingNet[i+1], 1, MPI_COMM_WORLD);
        	MPI_Recv(newPoints, elementsNumber, pointType, sortingNet[i+1], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	for(int i = 0, j = 0, k = 0; k < elementsNumber;) {
       			if ((*points)[i] < newPoints[j]) finalPoints[k++] = (*points)[i++];
       			else finalPoints[k++] = newPoints[j++];
        	}
        	*points = finalPoints;
    	}
    	else if (rank == sortingNet[i + 1]) {
        	Point *newPoints = new Point[elementsNumber], *finalPoints = new Point[elementsNumber];
        	MPI_Recv(newPoints, elementsNumber, pointType, sortingNet[i], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        	MPI_Send(*points, elementsNumber, pointType, sortingNet[i], 1, MPI_COMM_WORLD);
        	for (int i = elementsNumber - 1, j = elementsNumber - 1, k = elementsNumber - 1; k >= 0;) {
   				if (newPoints[j] < (*points)[i]) finalPoints[k--] = (*points)[i--];
   				else finalPoints[k--] = newPoints[j--];
        	}
        	*points = finalPoints;
    	}
	}
}

int main(int argc, char** argv) {
	if (argc < 3) {
    	cout << "n1 and n2 are needed" << endl;
    	return 0;
	}
	int n = stoi(argv[1]);
	int numprocs, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double t1, t2;
	int n1 = numprocs;
	int n2 = ceil(n / (float)n1);
	MPI_Barrier(MPI_COMM_WORLD);

	srand(time(NULL) + rank);

	createPointDatatype();

	Point *points = new Point[n2];
	generatePoints(points, n2, rank, n);
    
	vector<int> sortingNet = createSortingNet(0, numprocs);
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();

	distributedSort(&points, n2, rank, sortingNet);
	//sort(points, points + n2);
    
	MPI_Barrier(MPI_COMM_WORLD);
	t2 = MPI_Wtime();
	cout << "rank = " << rank << " sort time " << t2 - t1 << endl;
 
	MPI_Finalize();

	return 0;
}
