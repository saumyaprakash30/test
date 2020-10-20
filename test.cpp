#include<bits/stdc++.h>
// #include "mpi.h"
using namespace std;

void printFn(int a[],int n){
    for(int i=0;i<n;i++){
        cout<<a[i][0]<<endl;;
    }
}

int main(){
    int a[10][10];
    for(int i=0;i<10;i++)
       for(int j=0;j<10;j++)
            a[i][j] = i;
    
    printFn(&a,10);


}