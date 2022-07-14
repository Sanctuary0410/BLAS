#include<iostream>
#include <vector>
using namespace std;
#include "BLAS.h"
#include "myprint.h"

//≤‚ ‘dasum
/*int main()
{
	int len = 19;
	double arrA[19] = { 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1 };
	//double arrB[19] = { 2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2 };
	int inca = -1;
	double sum;
	//incb = 1;
	sum = dasum(19,arrA,inca);
	printArray(1, &sum);
	return 0;
}*/

/*≤‚ ‘daxpy*/
/*int main()
{
	int len = 19;
	double arrA[19] = { 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1 };
	double arrB[19] = { 2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2 };
	int inca = 1;
	int incb = 1;
	double sum;
	incb = 1;
	daxpy(19, 1, arrA, arrB, inca, incb);
	
	printArray(19, arrB);
	system("pause");
	return 0;
}*/

/*≤‚ ‘dcopy*/
/*int main()
{
	int len = 19;
	double arrA[19] = { 1,1,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1 };
	double arrB[19] = { 2,2,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2 };
	int inca = -2;
	int incb = -2;
	dcopy(10, arrA, arrB, inca, incb);
	printArray(19, arrB);
	system("pause");
	return 0;
}*/

/*≤‚ ‘ddot*/
/*int main()
{
	int len = 19;
	double arrA[19] = { 1,1,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1 };
	double arrB[19] = { 2,2,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2 };
	int inca = 1;
	int incb = 1;
	double sum;
	incb = 1;
	double dot=ddot(19, arrA, arrB, inca, incb);
	printArray(1, &dot);
	system("pause");
	return 0;
}*/

/*≤‚ ‘drot*/
/*int main()
{
	int len = 19;
	double arrA[19] = { 1,1,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1 };
	double arrB[19] = { 2,2,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2 };
	int inca = 1;
	int incb = 1;
	double sum;
	incb = 1;
	drot(19, arrA, arrB, 0, 1, inca, incb);
	printArray(19, arrB);
	system("pause");
	return 0;
}*/

/*≤‚ ‘drotg*/
/*int main()
{
	double x = 3, y = 4, c = 0, s = 0;
	drotg(&x, &y, &c, &s);
	printArray(1, &x);
	printArray(1, &y);
	printArray(1, &c);
	printArray(1, &s);
	system("pause");
	return 0;
}*/

/*c≤‚ ‘dscal*/
/*int main()
{
	int len = 19;
	double arrA[19] = { 1,1,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1 };
	//double arrB[19] = { 2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2 };
	int inca = 2;

	//incb = 1;
	dscal(19,2,arrA,inca);
	printArray(19, arrA);
	system("pause");
	return 0;
}*/

/*≤‚ ‘dswap*/
/*int main()
{
	int len = 19;
	double arrA[19] = { 1,1,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1 };
	double arrB[19] = { 2,2,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2 };
	int inca = -2;
	int incb = -2;
	dswap(10, arrA, arrB, inca, incb);
	printArray(19, arrB);
	system("pause");
	return 0;
}*/

/*c≤‚ ‘idamax*/
/*int main()
{
	int len = 19;
	double arrA[19] = { 1,10,2,0,1,0,1,0,1,0,1,0,1,0,8,0,1,0,1 };
	//double arrB[19] = { 2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2,0,2 };
	int inca = 2;
	double max;
	//incb = 1;
	max = idamax(10,arrA,inca);
	printArray(1, &max);
	system("pause");
	return 0;
}*/

/*c≤‚ ‘dnrm2*/
/*int main()
{
	double arrA[2] = {3,4};
	int inca = 1;
	double max;
	//incb = 1;
	double dnrom = dnrm2(2,arrA,inca);
	printArray(1, &dnrom);
	system("pause");
	return 0;
}*/

/*≤‚ ‘dger*/
/*int main()
{
	int len = 19, inca = 1, incb = 1, row = 5, col = 5;
	double arrA[5] = { 1,1,1,0,1 };
	double arrB[5] = { 2,2,2,0,2 };
	//double arrA[9] = { 1,1,1,0,1,0,1,0,1 };
	//double arrB[9] = { 2,2,2,0,2,0,2,0,2 };
	vector<vector<double> > a;
	a = { {0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0} };
	double alpha = 2;
	dger(arrA, arrB, alpha, a, row, col, inca, incb, 4);
	//er(double* x, double* y, double alpha, double** A, int m, int n, int incx, int incy)
	printMatrix(row, col, a);
	system("pause");
	return 0;
}*/

/*≤‚ ‘dgbmv*/
/*int main()
{
	int  inca = 2, incb =1, row = 4, col = 6, KL = 1, KU = 2;
	//double arrA[6] = { 1,1,1,0,1,0 };
	//double arrA[11] = { 1,0,1,0,1,0,0,0,1,0,0 };
	//double arrB[4] = { 0,0,0,0 };
	//double arrB[7] = { 0,0,0,0,0,0,0};
    //double arrA[4] = { 1,2,3,4 };
	double arrA[7] = { 1,0,2,0,3,0,4 };
	double arrB[6] = { 0,0,0,0,0,0};
	vector<vector<double> > a;
	//a = { {1,2,3,0,0,0},{1,2,3,4,0,0},{0,1,2,3,4,0},{0,0,1,2,3,4} };
	a= {  {0,0,1,2,3},
		  {0,1,2,3,4},
		  {0,1,2,3,4},
		  {0,1,2,3,4} };
	double alpha = 2,beta=1;
	char trans = 'C';
	dgbmv(a, arrA, arrB, alpha, beta, trans, row, col, KL, KU, 10, inca, incb);
	printArray(6, arrB);
	system("pause");
	return 0;
}*/

/*≤‚ ‘dtbsv*/
/*int main()
{
	int  n = 4, LDA = 10, incx = 2, k = 2;
	char trans = 'T', UPLO = 'L', DIAG = 'N';
	vector<vector<double> > a;
	
	//a = { {1,2,3,4},{1,2,3,4},{1,2,3,4},{1,2,3,4} };
    //a=U;
	//a = { {1,2,3},{2,3,4},{3,4,0},{4,0,0} };
	//a=L;
	a = { {0,0,1},
		  {0,1,2},
		  {1,2,3},
		  {2,3,4} };
	//double x[4] = { 1,2,3,4 };
	double x[7] = { 1,0,2,0,3,0,4 };
		  
	dtbsv(a, x, UPLO, trans, DIAG, n, k, LDA, incx);
	//(vector<vector<double>>& A, double* x, char UPLO, char trans, char diag, int n, int k, int LDA, int incx)
	printArray(7, x);
	system("pause");
	return 0;
}*/


/*≤‚ ‘dsyr2*/
/*int main()
{
	int len = 19, inca = 2, incb = 1, n = 5;
	double arrA[9] = { 1,0,0,0,1,0,0,0,1 };
	//double arrA[9] = { 1,0,1,0,1 };
	double arrB[5] = { 2,2,0,2,2 };
	vector<vector<double> > a;
	//a = { {2,2,0,0,0},{2,2,2,0,0},{0,2,2,2,0},{0,0,2,2,2},{0,0,0,2,2} };
	a = { {0,0,0,0,0},{0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0} };
	double alpha = 2;
	char UPLO = 'U';
	dsyr2(alpha, arrA, arrB, a, UPLO, n, inca, incb, 10);
	printMatrix(5, 5, a);
	system("pause");
	return 0;
}*/

/*≤‚ ‘dsymv*/
/*int main()
{
	int n=5, inca = -2, incb = -1;
	double alpha = 1, beta = 1;
	char UPLO = 'L';
	double arrA[9] = { 1,0,0,0,1,0,0,0,1 };
	//double arrA[5] = { 1,0,1,0,1 };
	//double arrB[5] = { 2,2,0,2,2 };
	double arrB[5] = { 1,0,0,0,1 };
	vector<vector<double> > a;
	a = { {2,2,0,0,0},{2,2,2,0,0},{0,2,2,2,0},{0,0,2,2,2},{0,0,0,2,2} };
	//a = { {0,0,0,0,0},{0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0} };
	
	dsymv(a, arrA, arrB, alpha, beta, UPLO, n, 10, inca, incb);
	printArray(5, arrB);
	system("pause");
	return 0;
}*/

/*≤‚ ‘dsymv*/
/*int main()
{
	int m = 4, n = 6,k=3,LDA=10,LDB=10,LDC=10;
	char transA = 'N', transB = 'T';
	double alpha = 1, beta = 1;
	vector<vector<double> > a,b,c;
	a = { {1,0,1},
		  {1,0,1}, 
		  {1,0,1},
		  {1,0,1} };
	//a = { {1, 1, 1, 1},{0 ,0, 0, 0},{1, 1, 1, 1} };
	//b = { {0,2,0,1,2,0},{0,2,0,1,2,0},{0,2,0,1,2,0} };
	b = { {0,0,0},
		  {2,2,2},
		  {0,0,0},
		  {1,1,1},
		  {2,2,2},
		  {0,0,0} };
	c = { {0,0,0,0,0,0},
		  {0,0,0,0,0,0} ,
		  {0,0,0,0,0,0} ,
		  {0,0,0,0,0,0} };
	dgemm(transA, transB, m, n, k, alpha,a,LDA,b,LDB,beta,c,LDC);
	printMatrix(4,6,c);
	system("pause");
	return 0;
}*/

/*≤‚ ‘dgemv*/
int main()
{
	int len = 19, inca = -2, incb = 1, row = 4, col = 6, KL = 1, KU = 2;
	//double arrA[6] = { 1,1,1,0,1,0 };
	//double arrB[4] = { 0,0,0,0 };
	double arrA[7] = { 0,0,1,0,1,0,1 };
	double arrB[6] = { 0,0,0,0,0,0 };
	vector<vector<double> > a;
	a = { {2,2,2,0,0,0},
		  {2,2,2,2,0,0},
		  {0,2,2,2,2,0},
		  {0,0,2,2,2,2} };
	double alpha = 2, beta = 1;
	char trans = 'T';
	dgemv(a, arrA, arrB, alpha, beta, trans, row, col, 10, inca, incb);
	printArray(6, arrB);
	return 0;
}

/*≤‚ ‘dtrmm*/
/*int main()
{
	int m = 4, n = 4,LDA=10,LDB=10;
	char SIDE='R',UPLO='L',transA = 'T',DIAG='N';
	double alpha = 1;
	vector<vector<double> > a,b;
	a = { {1,2,3,4},
		  {1,2,3,4},
		  {1,2,3,4},
		  {1,2,3,4} };
	b = { { 1,0,1,2}, 
		  { 1,0,1,2},
	      { 1,0,1,2},
	      { 1,0,1,2} };
	dtrmm(SIDE,UPLO,transA, DIAG, m, n, alpha,a,LDA,b,LDB);
	printMatrix(4,4,b);
	system("pause");
	return 0;
}*/

/*≤‚ ‘dtrsm*/
/*int main()
{
	int m = 4, n = 4, LDA = 10, LDB = 10;
	char SIDE = 'R', UPLO = 'L', transA = 'T', DIAG = 'U';
	double alpha = 1;
	vector<vector<double> > a, b;
	a = { {1,2,3,4},
		  {1,2,3,4},
		  {1,2,3,4},
		  {1,2,3,4} };
	b = { { 1,0,0,0},
		  { 0,1,0,0},
		  { 0,0,1,0},
		  { 0,0,0,1} };
	dtrsm(SIDE, UPLO, transA, DIAG, m, n, alpha, a, LDA, b, LDB);
	printMatrix(4, 4, b);
	system("pause");
	return 0;
}*/
