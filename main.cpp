#include<bits/stdc++.h>
#include "mpi.h"
using namespace std;

#define MASTER 0
#define PI 3.141592


struct complexNum{
	double real,imag;
	double value;
};

struct complexNum dft(int m,int n){
	double sum_r=0,sum_i=0;
	for(int i=0;i<ROW;i++){
		for (int j=0;j<COL;j++){
			double cosTerm = cos(2*PI*((double)((double)m*(double)i/(double)ROW) + (double)((double)n*(double)j/(double)COL)));
			double sinTerm = -1*sin(2*PI*((double)((double)m*(double)i/(double)ROW) + (double)((double)n*(double)j/(double)COL)));
			sum_r+=(double)img[i][j]*cosTerm;
			sum_i+=(double)img[i][j]*sinTerm;
		}
	}
	struct complexNum res;
	res.real = sum_r;
	res.imag = sum_i;
	res.value = sqrt(res.real*res.real+res.imag*res.imag)
	return res;


}
#define SIZE 200

int main(int argc, char *argv[]){
	int i,j;
	int total_proc;  
	int rank;        
	int n_per_proc,leftOut;
	int img[SIZE][SIZE];
	struct complexNum dftImg[SIZE][SIZE];

	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &total_proc);

	MPI_Comm_rank (MPI_COMM_WORLD,&rank);
	n_per_proc = SIZE/total_proc;
	leftOut = SIZE%total_proc;
	if(rank==MASTER){
		freopen("lena200.txt", "r", stdin);

		for(i=0;i<SIZE;i++){
			for(j=0;j<SIZE;j++){
				cin>>img[i][j];
			}
		}

	}
	MPI_Bcast(img, SIZE*SIZE, MPI_INT, 0, MPI_COMM_WORLD);
	struct complexNum temp[n_per_proc][SIZE];
	for(i=rank*n_per_proc;i<rank*n_per_proc+n_per_proc;i++){
		for(j=0;j<SIZE;j++){
			// struct complexNum temp1
			// temp1 = dft(i,j);
			// temp[i][j].real = temp1.real;
			// temp[i][j].imag = temp1.imag;
			temp[i-rank*n_per_proc][j] = dft(i,j);
		}
	}
	if(leftOut!=0 && rank==MASTER){
		for(i=n_per_proc*total_proc+1;i<SIZE;i++){
			for(j=0;j<SIZE;j++){
				
			}
		}
	}
	


	return 0;
}

