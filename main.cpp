#include <bits/stdc++.h>
#include "mpi.h"
using namespace std;

#define MASTER 0
#define PI 3.141592
#define SIZE 200
int img[SIZE][SIZE];

// struct complexNum{
//         double real,imag;
//         double value;
// };

// struct complexNum dft(int m,int n){
//         double sum_r=0,sum_i=0;
//         for(int i=0;i<SIZE;i++){
//                 for (int j=0;j<SIZE;j++){
//                         double cosTerm = cos(2*PI*((double)((double)m*(double)i/(double)SIZE) + (double)((double)n*(double)j/(double)SIZE)));
//                         double sinTerm = -1*sin(2*PI*((double)((double)m*(double)i/(double)SIZE) + (double)((double)n*(double)j/(double)SIZE)));
//                         sum_r+=(double)img[i][j]*cosTerm;
//                         sum_i+=(double)img[i][j]*sinTerm;
//                 }
//         }
//         struct complexNum res;
//         res.real = sum_r;
//         res.imag = sum_i;
//         res.value = sqrt(res.real*res.real+res.imag*res.imag);
//         return res;

// }

int main(int argc, char *argv[])
{
	int i, j, k, m, n;
	int total_proc;
	int rank;
	int n_per_proc, leftOut;
	int dftImg[SIZE][SIZE];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &total_proc);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	n_per_proc = SIZE / total_proc;
	leftOut = SIZE % total_proc;

	double start;

	if (rank == MASTER)
	{
		freopen("lena200.txt", "r", stdin);

		for (i = 0; i < SIZE; i++)
		{
			for (j = 0; j < SIZE; j++)
			{
				cin >> img[i][j];
			}
		}
	}
	// const int nitems=3;
	// int          blocklengths[3] = {1,1,1};
	// MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	// MPI_Datatype mpi_complexNum;
	// MPI_Aint     offsets[3];

	// offsets[0] = offsetof(complexNum, real);
	// offsets[1] = offsetof(complexNum, imag);
	// offsets[2] = offsetof(complexNum, value);

	// MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_complexNum);
	// MPI_Type_commit(&mpi_complexNum);

	start = MPI_Wtime();
	MPI_Status status;
	MPI_Bcast(img, SIZE * SIZE, MPI_INT, 0, MPI_COMM_WORLD);
	int temp[n_per_proc][SIZE];
	for (m = rank * n_per_proc; i < rank * n_per_proc + n_per_proc; i++)
	{
		for (n = 0; j < SIZE; j++)
		{
			double sum_r = 0, sum_i = 0;
			for (int i = 0; i < SIZE; i++)
			{
				for (int j = 0; j < SIZE; j++)
				{
					double cosTerm = cos(2 * PI * ((double)((double)m * (double)i / (double)SIZE) + (double)((double)n * (double)j / (double)SIZE)));
					double sinTerm = -1 * sin(2 * PI * ((double)((double)m * (double)i / (double)SIZE) + (double)((double)n * (double)j / (double)SIZE)));
					sum_r += (double)img[i][j] * cosTerm;
					sum_i += (double)img[i][j] * sinTerm;
				}
			}
			// struct complexNum res;
			// res.real = sum_r;
			// res.imag = sum_i;
			int value = sqrt(sum_r * sum_r + sum_i * sum_i);
			if (rank == MASTER)
			{
				// struct complexNum temp2;
				// temp2 = dft(i,j);
				dftImg[m][n] = value;
			}
			else
			{
				temp[m - rank * n_per_proc][n] = value;
			}
		}
	}
	if (rank != MASTER)
	{
		MPI_Send(&temp, n_per_proc * SIZE, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	if (rank == MASTER)
	{
		if (leftOut != 0)
		{
			int tempLeft[leftOut][SIZE];
			for (m = n_per_proc * total_proc; i < SIZE; i++)
			{
				for (n = 0; j < SIZE; j++)
				{
					double sum_r = 0, sum_i = 0;
					for (int i = 0; i < SIZE; i++)
					{
						for (int j = 0; j < SIZE; j++)
						{
							double cosTerm = cos(2 * PI * ((double)((double)m * (double)i / (double)SIZE) + (double)((double)n * (double)j / (double)SIZE)));
							double sinTerm = -1 * sin(2 * PI * ((double)((double)m * (double)i / (double)SIZE) + (double)((double)n * (double)j / (double)SIZE)));
							sum_r += (double)img[i][j] * cosTerm;
							sum_i += (double)img[i][j] * sinTerm;
						}
					}

					int value = sqrt(sum_r * sum_r + sum_i * sum_i);

					// tempLeft[i-(n_per_proc*total_proc)][j] = dft(i,j);
					// struct complexNum temp2
					// temp2 = dft(i,j);
					dftImg[m][n] = value;
				}
			}
		}
		for (int i = 1; i < total_proc; i++)
		{
			int tempF[n_per_proc][SIZE];
			MPI_Recv(&tempF, n_per_proc * SIZE, MPI_INT, i, 0, MPI_COMM_WORLD, status);
			for (j = i * n_per_proc; j < i * n_per_proc + n_per_proc; i++)
			{
				for (k = 0; k < SIZE; k++)
				{
					dftImg[j][k] = tempF[j - i * n_per_proc][k];
				}
			}
		}
		cout << MPI_Wtime() - start << endl;
	}

	return 0;
}