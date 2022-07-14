#include <stdio.h>
#include <vector>
using namespace std;

void printArray(double* arr,int len);
/*输出数组的前len个元素*/
void printArray(double* arr,int len)
{
        for (int i = 0; i < len; i++)
        {
                printf("%/n",arr[i]);
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
