#include<iostream>
#include <vector>
using namespace std;

void printArray(int len, double* arr);
void printMatrix(int m, int n, vector<vector<double> > A);


void printArray(int len, double* arr)
{
	for (int i = 0; i < len; i++)
	{
		cout << arr[i] << endl;
	}
}

void printMatrix(int m,int n, vector<vector<double>> A)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) {
			cout << A[i][j]<<',' ;
		}
		cout << endl;
	}
}
