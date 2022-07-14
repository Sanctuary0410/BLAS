#include<iostream>
#include <vector>
using namespace std;
#include <math.h>
#include "xerbal.h"
#define zero 0.00
#define one 1.00 


double dasum(int n, double* dx, int incx);

void daxpy(int n, double da, double* dx, double* dy, int incx, int incy);

void dcopy(int n, double* dx, double* dy, int incx, int incy);

double ddot(int n, double* dx, double* dy, int incx, int incy);

void drot(int n, double* dx, double* dy, double c, double s, int incx, int incy);

void drotg(double* da, double* db, double *c, double *s);

void dscal(int n, double da, double* dx, int  incx);

void dswap(int n, double* dx, double* dy, int incx, int incy);

double idamax(int n, double* dx, int incx);

double dnrm2(int n, double* x, int  incx);


void dger(double* x, double* y, double alpha, vector<vector<double> > &A, 
	int m, int n, int incx, int incy,int lda);

void dgemv(vector<vector<double> >& A, double* x, double* y, double alpha, double beta,
	char trans, int m, int n, int LDA, int incx, int incy);

void dgbmv(vector<vector<double> >& A, double* x, double* y, double alpha, double beta,
	char trans, int m, int n, int KL, int KU, int LDA, int incx, int incy);

void dtbsv(vector<vector<double>> &A, double* x, char UPLO, char trans, char diag, int n, int k, int LDA, int incx);

void dsyr2(double alpha, double* x, double* y, vector<vector<double>>& A,
	char UPLO, int n,int incx, int  incy,  int  LDA);

void dsymv(vector<vector<double> >& A, double* x, double* y, double alpha, double beta,
	char UPLO, int n, int LDA, int incx, int incy);

void dtbsv(vector<vector<double>>& A, double* x, char UPLO, char trans, char diag, int n, int k, int LDA, int incx);

void dgemm(char transA, char transB, int m, int n, int k, double alpha, vector < vector<double>>& A,
	int LDA, vector < vector<double>>& B, int LDB, double beta, vector < vector<double>>& C, int LDC);

void dtrmm(char SIDE, char UPLO, char transA, char DIAG, int m, int n, double alpha,
	vector < vector<double>>& A, int LDA, vector < vector<double>>& B, int LDB);

void dtrsm(char SIDE, char UPLO, char transA, char DIAG, int  m, int n, double alpha,
	vector<vector<double>>& A, int  LDA, vector<vector<double>>& B, int LDB);

//BLAS(L1)
//计算欧几里得1-范数
double dasum(int n, double* dx, int incx) {
	double dtemp = 0.0;
	int i, ix, m;
	if (n <= 0) {
		return 0;
	}
	if (incx == 1) {
		m = n % 6;
		if (m != 0) {
			for (i = 0; i < m; i++) {
				dtemp = dtemp + fabs(dx[i]);
			}
		}
		if (n >= 6) {
			for (i = m; i < n; i = i + 6)
				dtemp = dtemp + fabs(dx[i]) + fabs(dx[i + 1]) + fabs(dx[i + 2]) + fabs(dx[i + 3]) + fabs(dx[i + 4]) + fabs(dx[i + 5]);
		}
	}
	else {
		ix = 0;
		if (incx < 0) {
			ix = (-n + 1) * incx;
		}
		for (i = 0; i < n; i++) {
			dtemp = dtemp + fabs(dx[ix]);
			ix = ix + incx;
		}
	}
	return dtemp;
}

//将dx倍增加到dy上
void daxpy(int n, double da, double* dx, double* dy, int incx, int incy) {
	if ((n <= 0) || (da == 0)) {
		return;
	}
	int i;
	if ((incx == 1) & (incy == 1)) {
		int m = n % 4;
		if (m != 0) {
			for (i = 0; i < m; i++) {
				dy[i] = dy[i] + da * dx[i];
			}
		}
		if (n >= 4) {
			for (i = m; i < n; i = i + 4) {
				dy[i] = dy[i] + da * dx[i];
				dy[i + 1] = dy[i + 1] + da * dx[i + 1];
				dy[i + 2] = dy[i + 2] + da * dx[i + 2];
				dy[i + 3] = dy[i + 3] + da * dx[i + 3];
			}
		}
		return;
	}
	else {
		int ix = 0, iy = 0;
		if (incx < 0) {
			ix = (-n + 1) * incx;
		}
		if (incy < 0) {
			iy = (-n + 1) * incy;
		}
		for (i = 0; i < n; i++) {
			dy[iy] = dy[iy] + da * dx[ix];
			ix = ix + incx;
			iy = iy + incy;

		}
		return;
	}

}

//向量复制
void dcopy(int n, double* dx, double* dy, int incx, int incy) {
	if (n <= 0) {
		return;
	}
	int i;
	if ((incx == 0) & (incy == 0)) {
		int m = n % 7;
		if (m != 0) {
			for (i = 0; i < m; i++) {
				dy[i] = dx[i];
			}
		}
		if (n >= 7) {
			for (i = m; i < n; i = i + 7) {
				dy[i] = dx[i];
				dy[i + 1] = dx[i + 1];
				dy[i + 2] = dx[i + 2];
				dy[i + 3] = dx[i + 3];
				dy[i + 4] = dx[i + 4];
				dy[i + 5] = dx[i + 5];
				dy[i + 6] = dx[i + 6];
			}
		}
		return;
	}
	else {
		int ix = 0, iy = 0;
		if (incx < 0) {
			ix = (-n + 1) * incx;
		}
		if (incy < 0) {
			iy = (-n + 1) * incy;
		}
		for (i = 0; i < n; i++) {
			dy[iy] = dx[ix];
			ix = ix + incx;
			iy = iy + incy;
		}
	}
}

//求向量的点积
double ddot(int n, double* dx, double* dy, int incx, int incy) {
	double dtemp = 0;
	if (n <= 0) {
		return dtemp;
	}
	int i;
	if ((incx == 1) & (incy == 1)) {
		int m = n % 5;
		if (m != 0) {
			for (i = 0; i < m; i++) {
				dtemp = dtemp + dx[i] * dy[i];
			}
		}
		if (n >= 5) {
			for (i = m; i < n; i = i + 5) {
				dtemp = dtemp + dx[i] * dy[i] + dx[i + 1] * dy[i + 1] + dx[i + 2] * dy[i + 2] + dx[i + 3] * dy[i + 3] + dx[i + 4] * dy[i + 4];
			}
		}
	}
	else {
		int ix = 0, iy = 0;
		if (incx < 0) {
			ix = (-n + 1) * incx;
		}
		if (incy < 0) {
			iy = (-n + 1) * incy;
		}
		for (i = 0; i < n; i++) {
			dtemp = dtemp + dx[ix] * dy[iy];
			ix = ix + incx;
			iy = iy + incy;
		}

	}
	return dtemp;
}

//实现平面交换
void drot(int n, double *dx, double *dy, double c, double s, int incx, int incy) {
	double dtemp;
	if (n <= 0) {
		return;
	}
	int i;
	if ((incx == 1) & (incy == 1)) {
		for (i = 0; i < n; i++) {
			dtemp = c * dx[i] + s * dy[i];
			dy[i] = c * dy[i] - s * dx[i];
			dx[i] = dtemp;
		}
		return;
	}
	else {
		int ix = 0, iy = 0;
		if (incx < 0) {
			ix = (-n + 1) * incx;
		}
		if (incy < 0) {
			iy = (-n + 1) * incy;
		}
		for (i = 0; i < n; i++) {
			dtemp = c * dx[ix] + s * dy[iy];
			dy[iy] = c * dy[iy] - s * dx[ix];
			dx[ix] = dtemp;
			ix = ix + incx;
			iy = iy + incy;
		}
		return;
	}
}

//建立平面交换
void drotg(double* da, double* db, double* c, double* s) {
	double roe = *db, r, z;
	if (fabs(*da) > fabs(*db)) {
		roe = *da;
	}
	double scale = fabs(*da) + fabs(*db);
	if (scale != 0) {
		r = scale * sqrt(pow(*da / scale, 2) + pow(*db / scale, 2));
		if (roe > 0) {
			r = r;
		}
		else {
			r = -r;
		}
		*c = *da / r;
		*s = *db / r;
		z = 1;
		if (fabs(*da) > fabs(*db)) {
			z = *s;
		}
		else if (c != 0) {
			z = 1.0 / *c;
		}
	}
	else {
		*c = 1.0;
		*s = 0.0;
		r = 0.0;
		z = 0.0;
	}
	*da = r;
	*db = z;
}

//向量缩放
void dscal(int n, double da, double* dx, int  incx) {
	if (n <= 0) {
		return;
	}
	int i;
	if (incx == 1) {
		int m = n % 5;
		if (m != 0) {
			for (i = 0; i < m; i++) {
				dx[i] = da * dx[i];
			}
		}
		if (n >= 5) {
			for (i = m; i < n; i = i + 5) {
				dx[i] = da * dx[i];
				dx[i + 1] = da * dx[i + 1];
				dx[i + 2] = da * dx[i + 2];
				dx[i + 3] = da * dx[i + 3];
				dx[i + 4] = da * dx[i + 4];
			}
		}
	}
	else {
		int ix = 0;
		if (incx < 0) {
			ix = (-n + 1) * incx;
		}
		for (i = 0; i < n; i++) {
			dx[ix] = da * dx[ix];
			ix = ix + incx;
		}
	}
	return;
}

//向量交换
void dswap(int n, double* dx, double* dy, int incx, int incy) {
	if (n <= 0) {
		return;
	}
	double dtemp;
	int i;
	if ((incx == 1) & (incy == 1)) {
		int m = n % 3;
		if (m != 0) {
			for (i = 0; i < m; i++) {
				dtemp = dx[i];
				dx[i] = dy[i];
				dy[i] = dtemp;
			}
		}
		if (n >= 3) {
			for (i = m; i < n; i = i + 3) {
				dtemp = dx[i];
				dx[i] = dy[i];
				dy[i] = dtemp;
				dtemp = dx[i + 1];
				dx[i + 1] = dy[i + 1];
				dy[i + 1] = dtemp;
				dtemp = dx[i + 2];
				dx[i + 2] = dy[i + 2];
				dy[i + 2] = dtemp;
			}
		}
	}
	else {
		int ix = 0, iy = 0;
		if (incx < 0) {
			ix = (-n + 1) * incx;
		}
		if (incy < 0) {
			iy = (-n + 1) * incy;
		}
		for (i = 0; i < n; i++) {
			dtemp = dx[ix];
			dx[ix] = dy[iy];
			dy[iy] = dtemp;
			ix = ix + incx;
			iy = iy + incy;
		}
	}
	return;
}

//求欧几里得∞-范数
double idamax(int n, double* dx, int incx) {
	int idamax = 0, i;
	if (n < 1) {
		return 0;
	}
	double dmax = fabs(dx[0]);
	if (n == 1) {
		return dmax;
	}
	if (incx == 1) {
		for (i = 1; i < n; i++) {
			if (fabs(dx[i]) > dmax) {
				idamax = i;
				dmax = fabs(dx[i]);
			}
		}
	}
	else {
		int ix = 0;
		if (incx < 0) {
			ix = (-n + 1) * incx;
		}
		dmax = fabs(dx[ix]);
		ix = ix + incx;
		for (i = 1; i < n; i++) {
			if (fabs(dx[ix]) > dmax) {
				idamax = i;
				dmax = fabs(dx[ix]);
			}
			ix = ix + incx;
		}
	}
	return dmax;
}

//求欧几里得2-范数
double dnrm2(int n, double* x, int  incx) {
	double dnrm2 = 0;
	if (n <= 0) return dnrm2;
	int i = 0, ix, j, next = 30;
	double sum = 0.00, xmax, temp, cutlo = 0.00001, cuthi = 10000, hitest;

    /*cutlo = maximum of  dsqrt(u / eps)  over all known machines.
	cuthi = minimum of  dsqrt(v)      over all known machines.
    where
	eps = smallest no.such that eps + 1..gt. 1.
	u = smallest positive no.   (underflow limit)
	v = largest  no.            (overflow  limit)*/
	if (incx < 0) i = (-n + 1) * incx;
	ix = 1;
	while (ix <= n) {
		if (next == 30) {
			if (fabs(x[i]) > cutlo) {
				hitest = cuthi / n;
				while (ix <= n) {
					if (fabs(x[i]) < hitest) {
						sum = sum + x[i] * x[i];
						i = i + incx;
						ix = ix + 1;
					}
					else break;
					if (ix > n) return dnrm2 = sqrt(sum);
				}
				next = 110;
				sum = sum / (x[i] * x[i]);
				xmax = x[i];
				temp = x[i] / xmax;
				sum = sum + temp * temp;
				ix = ix + 1;
				i = i + incx;
			}
			else {
				next = 50;
				xmax = 0;
				if (x[i] == zero) {
					ix = ix + 1;
					i = i + incx;
				}
				else {
					next = 70;
					xmax = x[i];
					temp = x[i] / xmax;
					sum = sum + temp * temp;
					ix = ix + 1;
					i = i + incx;
				}
			}
		}
		else if (next == 50) {
			if (x[i] == zero) {
				ix = ix + 1;
				i = i + incx;
			}
			else if (fabs(x[i] > cutlo)) {
				hitest = cuthi / n;
				//j = ix;
				while (ix <= n) {
					if (fabs(x[i]) < hitest) {
						sum = sum + x[i] * x[i];
						i = i + incx;
						ix = ix + 1;
					}
					else break;
					if (ix > n) return dnrm2 = sqrt(sum);
				}
				next = 110;
				sum = sum / (x[i] * x[i]);
				xmax = x[i];
				temp = x[i] / xmax;
				sum = sum + temp * temp;
				ix = ix + 1;
				i = i + incx;
			}
			else {
				xmax = x[i];
				temp = x[i] / xmax;
				sum = sum + temp * temp;
				ix = ix + 1;
				i = i + incx;
			}
		}
		else if (next == 70) {
			if (fabs(x[i] > cutlo)) {
				sum = sum * xmax * xmax;
				hitest = cuthi / n;
				//j = ix;
				while (ix <= n) {
					if (fabs(x[i]) < hitest) {
						sum = sum + x[i] * x[i];
						i = i + incx;
						ix = ix + 1;
					}
					else break;
					if (ix > n) return dnrm2 = sqrt(sum);
				}
				next = 110;
				sum = sum / (x[i] * x[i]);
				xmax = x[i];
				temp = x[i] / xmax;
				sum = sum + temp * temp;
				ix = ix + 1;
				i = i + incx;
			}
			else {
				if (x[i] <= xmax) {
					temp = x[i] / xmax;
					sum = sum + temp * temp;
					ix = ix + 1;
					i = i + incx;
				}
				else {
					temp = x[i] / xmax;
					sum = 1 + sum * temp * temp;
					ix = ix + 1;
					i = i + incx;
				}
			}
		}
		else if (next == 110) {
			if (x[i] <= xmax) {
				temp = x[i] / xmax;
				sum = sum + temp * temp;
				ix = ix + 1;
				i = i + incx;
			}
			else {
				temp = x[i] / xmax;
				sum = 1 + sum * temp * temp;
				ix = ix + 1;
				i = i + incx;
			}

		}
	}
	return dnrm2 = xmax * sqrt(sum);

}

//BLAS(L2)
/*普通矩阵的单位校正A=alpha*xy'+A*/
void dger(double* x, double* y, double alpha, vector<vector<double> >& A, int m, int n, int incx, int incy, int LDA) {
	int info = 0;
	if (m < 0) info = 6;

	else if (n < 0) info = 7;

	else if (incx == 0) info = 8;

	else if (incy == 0) info = 9;

	else if (LDA < 0) info = 10;

	else if (info != 0) {
		char name[] = "dger";
		xerbla(name, info);
		return;
	}

	if ((m == 0) || (n == 0) || (alpha == 0)) return;

	int ix, j, i;
	double temp;
	if (incx > 0) {
		ix = 0;
	}
	else {
		ix = (1 - n) * incx;
	}
	if (incy == 1) {
		for (i = 0; i < m; i++) {
			if (x[ix] != 0) {
				temp = alpha * x[ix];
				for (j = 0; j < n; j++) {
					A[i][j] = A[i][j] + y[j] * temp;
				}
			}
			ix = ix + incx;
		}
	}
	else {
		int jy, ky;
		if (incy > 0) {
			ky = 0;
		}
		else {
			ky = (1 - m) * incy;
		}
		for (i = 0; i < m; i++) {
			if (x[ix] != 0) {
				temp = alpha * x[ix];
				jy = ky;
				for (j = 0; j < n; j++) {
					A[i][j] = A[i][j] + y[jy] * temp;
					jy = jy + incy;
				}
			}
			ix = ix + incx;
		}
	}
	return;
}

/*普通带状矩阵与向量乘积*/
void dgbmv(vector<vector<double> >& A, double* x, double* y, double alpha, double beta,
	char trans, int m, int n, int KL, int KU, int LDA, int incx, int incy) {
	/* TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.

	   TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.

	   TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.*/
	   /*The following program segment will transfer a band matrix
	   from conventional full matrix storage to band storage :
		   for (i = 0; i < n; i++) {
			   m = ku - i;
			   for (j = i; j <= min(n, j + k); j++) {
				   A[i][j+m]=matrix[i][j];
			   }
		   }使对角线元素在同一列
	   */
	int info = 0;
	if ((trans != 'N') & (trans != 'T') & (trans != 'C'))info = 6;

	else if (m < 0) info = 7;

	else if (n < 0) info = 8;

	else if (KL < 0)info = 9;

	else if (KU < 0)info = 10;

	else if (LDA < (KL + KU + 1))info = 11;

	else if (incx == 0) info = 12;

	else if (incy == 0) info = 13;

	if (info != 0) {
		char name[] = "dgbmv";
		xerbla(name, info);
		return;
	}
	if ((m == 0) || (n == 0) || ((alpha == zero) & (beta == zero))) {
		return;
	}
	int lenx, leny, kx, ky, i, j;
	double temp;
	if (trans == 'N') {
		lenx = n;
		leny = m;
	}
	else {
		lenx = m;
		leny = n;
	}
	if (incx > 0) kx = 0;
	else kx = (1 - lenx) * incx;
	if (incy > 0) ky = 0;
	else ky = (1 - leny) * incy;

	/*First form y:=beta*y*/
	if (beta != one) {
		if (incy == 1) {
			if (beta == zero) {
				for (i = 0; i < leny; i++) y[i] = 0;
			}
			else {
				for (i = 0; i < leny; i++) y[i] = beta * y[i];
			}
		}
		else {
			int iy = ky;
			if (beta == zero) {
				for (i = 0; i < leny; i++) {
					y[iy] = 0;
					iy = iy + incy;
				}
			}
			else {
				for (i = 0; i < leny; i++) {
					y[iy] = beta * y[iy];
					iy = iy + incy;
				}
			}
		}
	}

	if (alpha == zero)return;
	if (trans == 'N') {
		/*Form  y := alpha*A*x + y.*/
		int iy = ky;
		if (incx == 1) {
			for (i = 0; i < m; i++) {
				temp = zero;
				int l = KU - i;
				for (j = MAX(0, i - KL); j <= MIN(n - 1, i + KU); j++) {
					temp = temp + A[i][j + l] * x[j];
				}
				y[iy] = y[iy] + alpha * temp;
				iy = iy + incy;
			}
		}
		else {
			int jy = ky;
			for (i = 0; i < m; i++) {
				temp = zero;
				int max = MAX(0, i - KL);
				int jx = kx + incx * max;
				int l = KU - i;
				for (j = max; j <= MIN(n - 1, i + KU); j++) {
					temp = temp + A[i][j + l] * x[jx];
					jx = jx + incx;
				}
				y[jy] = y[jy] + alpha * temp;
				jy = jy + incy;
			}
		}

	}
	else {
		/*Form  y := alpha*A'*x + y.*/
		if (incx == 1) {
			for (i = 0; i < m; i++) {
				int l = KU - i;

				if (x[i] != zero) {

					temp = alpha * x[i];
					int max = MAX(0, i - KL);
					int jy = ky + incy * max;
					for (j = max; j < MIN(n, i + KU + 1); j++) {
						y[jy] = y[jy] + temp * A[i][j + l];
						jy = jy + incy;
					}
				}
			}
		}
		else {
			int ix = kx;
			for (i = 0; i < m; i++) {
				int l = KU - i;
				if (x[ix] != zero) {
					temp = alpha * x[ix];
					int max = MAX(0, i - KL);
					int jy = ky + incy * max;
					for (j = max; j < MIN(n, i + KU + 1); j++) {
						y[jy] = y[jy] + temp * A[i][j + l];
						jy = jy + incy;
					}
				}
				ix = ix + incx;
			}

		}

	}
	return;
}

/*普通矩阵与向量乘积y := alpha*A*x + beta*y*/
void dgemv(vector<vector<double> >& A, double* x, double* y, double alpha, double beta,
	char trans, int m, int n, int LDA, int incx, int incy) {
	/* TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.

	   TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.

	   TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.*/
	int info = 0;
	if ((trans != 'N') & (trans != 'T') & (trans != 'C'))info = 1;

	else if (m < 0) info = 2;

	else if (n < 0) info = 3;

	else if (LDA < MAX(1, m))info = 6;

	else if (incx == 0) info = 8;

	else if (incy == 0) info = 11;

	if (info != 0) {
		char name[] = "dgemv";
		xerbla(name, info);
		return;
	}
	if ((m == 0) || (n == 0) || ((alpha == zero) & (beta == zero))) {
		return;
	}
	int lenx, leny, kx, ky, i, j;
	double temp;
	if (trans == 'N') {
		lenx = n;
		leny = m;
	}
	else {
		lenx = m;
		leny = n;
	}
	if (incx > 0) kx = 0;
	else kx = (1 - lenx) * incx;
	if (incy > 0) ky = 0;
	else ky = (1 - leny) * incy;

	/*First form y:=beta*y*/
	if (beta != one) {
		if (incy == 1) {
			if (beta == zero) {
				for (i = 0; i < leny; i++) y[i] = 0;
			}
			else {
				for (i = 0; i < leny; i++) y[i] = beta * y[i];
			}
		}
		else {
			int iy = ky;
			if (beta == zero) {
				for (i = 0; i < leny; i++) {
					y[iy] = 0;
					iy = iy + incy;
				}
			}
			else {
				for (i = 0; i < leny; i++) {
					y[iy] = beta * y[iy];
					iy = iy + incy;
				}
			}
		}
	}

	if (alpha == zero)return;
	if (trans == 'N') {
		/*Form  y := alpha*A*x + y.*/
		int iy = ky;
		if (incx == 1) {
			for (i = 0; i < m; i++) {
				temp = zero;
				for (j = 0; j < n; j++) {
					temp = temp + A[i][j] * x[j];
				}
				y[iy] = y[iy] + alpha * temp;
				iy = iy + incy;
			}
		}
		else {
			int jy = ky;
			for (i = 0; i < m; i++) {
				temp = zero;
				int jx = kx;
				for (j = 0; j < n; j++) {
					temp = temp + A[i][j] * x[jx];
					jx = jx + incx;
				}
				y[jy] = y[jy] + alpha * temp;
				jy = jy + incy;
			}
		}

	}
	else {
		/*Form  y := alpha*A'*x + y.*/
		int jx = kx;
		if (incy == 1) {
			for (i = 0; i < m; i++) {
				if (x[jx] != zero) {
					temp = alpha * x[jx];
					for (j = 0; j < n; j++) {
						y[j] = y[j] + temp * A[i][j];
					}
				}
				jx = jx + incx;
			}
		}
		else {
			for (i = 0; i < m; i++) {
				if (x[jx] != zero) {
					int jy = ky;
					temp = alpha * x[jx];
					for (j = 0; j < n; j++) {
						y[jy] = y[jy] + temp * A[i][j];
						jy = jy + incy;
					}
				}
				jx = jx + incx;
			}
		}

	}
	return;
}

/*实数对称矩阵的2秩校正A:=alpha*xy'+alpha*yx'+A*/
void dsyr2(double alpha, double* x, double* y, vector<vector<double>>& A,
	char UPLO, int n, int incx, int  incy, int  LDA) {
	int info = 0;
	if ((UPLO != 'U') & (UPLO != 'L')) info = 5;
	else if (n < 0)info = 6;
	else if (incx == 0)info = 7;
	else if (incy == 0)info = 8;
	else if (LDA < MAX(1, n))info = 9;
	if (info != 0) {
		char name[] = "dsyr2";
		xerbla(name, info);
		return;
	}
	if ((n == 0) || (alpha == 0)) return;
	int kx, ky, jx, jy, i, j;
	double temp1, temp2;
	if ((incx != 1) || (incy != 1)) {
		if (incx > 0) kx = 0;
		else kx = (1 - n) * incx;
		if (incy > 0) ky = 0;
		else ky = (1 - n) * incy;
		jx = kx;
		jy = ky;
	}
	if (UPLO == 'U') {
		if ((incx == 1) & (incy == 1)) {
			for (i = 0; i < n; i++) {
				if ((x[i] != 0) || (y[i] != 0)) {
					temp1 = alpha * x[i];
					temp2 = alpha * y[i];
					for (j = i; j < n; j++) {
						A[i][j] = A[i][j] + temp1 * y[j] + temp2 * x[j];
					}
				}
			}
		}
		else {
			int ix, iy;

			for (i = 0; i < n; i++) {
				if ((x[jx] != 0) || (y[jy] != 0)) {
					temp1 = alpha * x[jx];
					temp2 = alpha * y[jy];
					ix = jx;
					iy = jy;
					for (j = i; j < n; j++) {
						A[i][j] = A[i][j] + temp1 * y[iy] + temp2 * x[ix];
						ix = ix + incx;
						iy = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
			}
		}
	}
	else {
		if ((incx == 1) & (incy == 1)) {
			for (i = 0; i < n; i++) {
				if ((x[i] != 0) || (y[i] != 0)) {
					temp1 = alpha * x[i];
					temp2 = alpha * y[i];
					for (j = 0; j < i + 1; j++) {
						A[i][j] = A[i][j] + temp1 * y[j] + temp2 * x[j];
					}
				}
			}
		}
		else {
			int ix, iy;
			for (i = 0; i < n; i++) {
				if ((x[jx] != 0) || (y[jy] != 0)) {
					temp1 = alpha * x[jx];
					temp2 = alpha * y[jy];
					ix = kx;
					iy = ky;
					for (j = 0; j < i + 1; j++) {
						A[i][j] = A[i][j] + temp1 * y[iy] + temp2 * x[ix];
						ix = ix + incx;
						iy = iy + incy;
					}
				}
				jx = jx + incx;
				jy = jy + incy;
			}
		}
	}
	return;
}

/*对称矩阵与向量乘积y := alpha*A*x + beta*y*/
void dsymv(vector<vector<double> >& A, double* x, double* y, double alpha, double beta,
	char UPLO, int n, int LDA, int incx, int incy) {
	int info = 0;
	if ((UPLO != 'U') & (UPLO != 'L'))info = 6;

	else if (n < 0) info = 7;

	else if (LDA < MAX(1, n))info = 8;

	else if (incx == 0) info = 9;

	else if (incy == 0) info = 10;

	if (info != 0) {
		char name[] = "dsymv";
		xerbla(name, info);
		return;
	}
	if ((n == 0) || ((alpha == zero) & (beta == zero))) {
		return;
	}

	int kx, ky, i, j;
	if (incx > 0) kx = 0;
	else kx = (1 - n) * incx;
	if (incy > 0) ky = 0;
	else ky = (1 - n) * incy;

	/*First form y:=beta*y*/
	if (beta != one) {
		if (incy == 1) {
			if (beta == zero) {
				for (i = 0; i < n; i++) y[i] = zero;
			}
			else {
				for (i = 0; i < n; i++) y[i] = beta * y[i];
			}
		}
		else {
			int iy = ky;
			if (beta == zero) {
				for (i = 0; i < n; i++) {
					y[iy] = 0;
					iy = iy + incy;
				}
			}
			else {
				for (i = 0; i < n; i++) {
					y[iy] = beta * y[iy];
					iy = iy + incy;
				}
			}
		}
	}

	if (alpha == zero)return;
	double temp1, temp2;
	if (UPLO == 'U') {
		if ((incx == 1) & (incy == 1)) {
			for (i = 0; i < n; i++) {
				temp1 = alpha * x[i];
				temp2 = zero;
				y[i] = y[i] + temp1 * A[i][i];
				for (j = i + 1; j < n; j++) {
					y[j] = y[j] + temp1 * A[i][j];
					temp2 = temp2 + A[i][j] * x[j];
				}
				y[i] = y[i] + alpha * temp2;
			}
		}
		else {
			int jx = kx;
			int jy = ky;
			int ix, iy;
			for (i = 0; i < n; i++) {
				temp1 = alpha * x[jx];
				temp2 = zero;
				y[jy] = y[jy] + temp1 * A[i][i];
				ix = jx + incx;
				iy = jy + incy;
				for (j = i + 1; j < n; j++) {
					y[iy] = y[iy] + temp1 * A[i][j];
					temp2 = temp2 + A[i][j] * x[ix];
					ix = ix + incx;
					iy = iy + incy;
				}
				y[jy] = y[jy] + alpha * temp2;
				jx = jx + incx;
				jy = jy + incy;
			}
		}
	}

	else {
		if ((incx == 1) & (incy == 1)) {
			for (i = 0; i < n; i++) {
				temp1 = alpha * x[i];
				temp2 = zero;
				for (j = 0; j < i; j++) {
					y[j] = y[j] + temp1 * A[i][j];
					temp2 = temp2 + A[i][j] * x[j];
				}
				y[i] = y[i] + temp1 * A[i][i] + alpha * temp2;
			}
		}
		else {
			int jx = kx;
			int ix, iy;
			for (i = 0; i < n; i++) {
				iy = ky;
				ix = kx;
				temp1 = alpha * x[jx];
				temp2 = zero;
				for (j = 0; j < i; j++) {
					y[iy] = y[iy] + temp1 * A[i][j];
					temp2 = A[i][j] * x[ix];
					iy = iy + incy;
					ix = ix + incx;
				}
				y[iy] = y[iy] + temp1 * A[i][i] + temp2 * alpha;
				jx = jx + incx;
			}
		}
	}
	return;
}

/*三角带状矩阵求逆后与向量乘积*/
void dtbsv(vector<vector<double>>& A, double* x, char UPLO, char trans, char diag, int n, int k, int LDA, int incx) {
	/*The following program segment will transfer an upper
	triangular band matrix from conventional full matrix storage
	to band storage :
		for (i = 0; i < n; i++) {
			m = - i;
			for (j = i; j <= min(n, j + k); j++) {
				A[i][j+m]=matrix[i][j];
			}
		}对角线元素存储在第一列

	The following program segment will transfer an lower
	triangular band matrix from conventional full matrix storage
	to band storage :
		for (i = 0; i < n; i++) {
			m = k - i;
			for (j = i; j <= min(n, j + k); j++) {
				A[i][j+m]=matrix[i][j];
			}
		}对角形元素存储在最后一列
	 */

	int info = 0;
	if ((UPLO != 'U') & (UPLO != 'L')) info = 3;
	else if ((trans != 'N') & (trans != 'T') & (trans != 'C'))info = 4;
	else if ((diag != 'U') & (diag != 'N')) info = 5;
	else if (n < 0)info = 6;
	else if (k < 0)info = 7;
	else if (LDA < k + 1)info = 8;
	else if (incx == 0)info = 9;

	if (info != 0) {
		char name[] = "dtbsv";
		xerbla(name, info);
		return;
	}

	if (n == 0)return;

	bool NOUNIT = (diag == 'N');

	int kx;
	if (incx < 0)kx = (1 - n) * incx;
	else if (incx != 1) kx = 0;

	int i, j, l;
	if (trans == 'N') {
		//Form  x : = inv(A) * x.
		if (UPLO == 'U') {
			if (incx == 1) {
				for (i = n - 1; i >= 0; i = i - 1) {
					l = -i;
					for (j = MIN(i + k, n - 1); j > i; j = j - 1) {
						x[i] = x[i] - x[j] * A[i][j + l];
					}
					if (NOUNIT) x[i] = x[i] / A[i][i + l];
				}
			}
			else {
				kx = kx + (n - 1) * incx;
				int ix = kx;
				for (i = n - 1; i >= 0; i = i - 1) {
					l = -i;
					int min = MIN(i + k, n - 1);
					int jx = kx - incx * (n - 1 - min);
					for (j = min; j > i; j = j - 1) {
						x[ix] = x[ix] - x[jx] * A[i][j + l];
						jx = jx - incx;
					}
					if (NOUNIT) x[ix] = x[ix] / A[i][i + l];
					ix = ix - incx;
				}
			}
		}
		else {
			if (incx == 1) {
				for (i = 0; i < n; i++) {
					l = k - i;
					for (j = MAX(0, i - k); j < i; j++) {
						x[i] = x[i] - x[j] * A[i][j + l];
					}
					if (NOUNIT) x[i] = x[i] / A[i][i + l];
				}
			}
			else {
				int ix = kx;
				for (i = 0; i < n; i++) {
					l = k - i;
					int max = MAX(0, i - k);
					int jx = kx + incx * max;
					for (j = max; j < i; j++) {
						x[ix] = x[ix] - x[jx] * A[i][j + l];
						jx = jx + incx;
					}
					if (NOUNIT) x[ix] = x[ix] / A[i][i + l];
					ix = ix + incx;
				}
			}
		}
	}
	else {
		//Form  x : = inv(A') * x.
		if (UPLO == 'U') {
			if (incx == 1) {
				for (i = 0; i < n; i++) {
					int l = -i;
					if (NOUNIT) x[i] = x[i] / A[i][i + l];
					for (j = i + 1; j <= MIN(n - 1, i + k); j++) {
						x[j] = x[j] - x[i] * A[i][j + l];
					}
				}
			}
			else {
				int ix = kx;
				for (i = 0; i < n; i++) {
					int l = -i;
					if (NOUNIT) x[ix] = x[ix] / A[i][i + l];
					int jx = ix + incx;
					for (j = i + 1; j <= MIN(n - 1, i + k); j++) {
						x[jx] = x[jx] - x[ix] * A[i][j + l];
						jx = jx + incx;
					}
					ix = ix + incx;
				}
			}
		}
		else {
			if (incx == 1) {
				for (i = n - 1; i >= 0; i = i - 1) {
					int l = k - i;
					if (NOUNIT) x[i] = x[i] / A[i][i + l];
					for (j = i - 1; j >= MAX(0, i - k); j = j - 1) {
						x[j] = x[j] - x[i] * A[i][j + l];
					}
				}
			}
			else {
				int ix = kx + incx * (n - 1);
				for (i = n - 1; i >= 0; i = i - 1) {
					int l = k - i;
					int jx = ix - incx;
					if (NOUNIT) x[ix] = x[ix] / A[i][i + l];
					for (j = i - 1; j >= MAX(0, i - k); j = j - 1) {
						x[jx] = x[jx] - x[ix] * A[i][j + l];
						jx = jx - incx;
					}
					ix = ix - incx;
				}
			}
		}

	}
}

//BLAS(L3)
/*普通矩阵乘法*/
void dgemm(char transA, char transB, int m, int n, int k, double alpha, vector < vector<double>>& A,
	int LDA, vector < vector<double>>& B, int LDB, double beta, vector < vector<double>>& C, int LDC) {
	/*Form C := alpha*op( A )*op( B ) + beta*C*/
	bool notA = (transA == 'N');
	bool notB = (transB == 'N');
	int nrowA = m, ncolA = k, nrowB = k;
	if (!notA) {
		nrowA = k;
		ncolA = m;
	}

	if (!notB)nrowB = n;

	/*Test the input parameters.*/
	int info = 0;
	if ((!notA) & (transA != 'C') & (transA != 'T'))info = 1;
	else if ((!notB) & (transB != 'C') & (transB != 'T'))info = 2;
	else if (m < 0)info = 3;
	else if (n < 0)info = 4;
	else if (k < 0)info = 5;
	else if (LDA < MAX(1, nrowA))info = 8;
	else if (LDB < MAX(1, nrowB))info = 10;
	else if (LDC < MAX(1, m))info = 13;

	if (info != 0) {
		char name[] = "dgemm";
		xerbla(name, info);
		return;
	}
	int i, j, l;
	if ((m == 0) || (n == 0) || (((alpha == zero) || (k == 0)) & (beta == 0)))return;
	if (alpha == zero) {
		if (beta == zero) {
			for (i = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					C[i][j] = zero;
				}
			}
		}
		else {
			for (i = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					C[i][j] = beta * C[i][j];
				}
			}
		}
		return;
	}

	double temp;
	if (notB) {
		for (i = 0; i < m; i++) {
			if (beta == zero) {
				for (j = 0; j < n; j++) {
					C[i][j] = zero;
				}
			}
			else if (beta != one) {
				for (j = 0; j < n; j++) {
					C[i][j] = beta * C[i][j];
				}
			}
		}
		if (notA) {
			/* Form  C := alpha*A*B + beta*C.*/
			for (i = 0; i < m; i++) {
				for (l = 0; l < k; l++) {
					if (A[i][l] != 0) {
						temp = alpha * A[i][l];
						for (j = 0; j < n; j++) {
							C[i][j] = C[i][j] + temp * B[l][j];
						}
					}
				}
			}
			return;
		}
		else {
			/*Form  C := alpha*A'*B + beta*C*/
			for (l = 0; l < k; l++) {
				for (i = 0; i < m; i++) {
					if (A[l][i] != 0) {
						temp = alpha * A[l][i];
						for (j = 0; j < n; j++) {
							C[i][j] = C[i][j] + temp * B[l][j];
						}
					}

				}

			}
			return;
		}
	}
	else {
		if (notA) {
			/*Form  C := alpha*A*B' + beta*C*/
			for (i = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					temp = zero;
					for (l = 0; l < k; l++) {
						temp = temp + A[i][l] * B[j][l];
					}
					if (beta == zero) {
						C[i][j] = temp * alpha;
					}
					else {
						C[i][j] = temp * alpha + beta * C[i][j];
					}
				}
			}
		}
		else {
			/*Form  C := alpha*A'*B' + beta*C*/
			for (i = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					temp = zero;
					for (l = 0; l < k; l++) {
						temp = temp + A[l][i] * B[j][l];
					}
					if (beta == zero) {
						C[i][j] = temp * alpha;
					}
					else {
						C[i][j] = temp * alpha + beta * C[i][j];
					}

				}

			}
		}
	}
	return;
}

/*三角矩阵与矩阵乘法*/
void dtrmm(char SIDE, char UPLO, char transA, char DIAG, int m, int n, double alpha,
	vector < vector<double>>& A, int LDA, vector < vector<double>>& B, int LDB) {
	/*Form B := alpha*op( A )*B,   or   B := alpha*B*op( A )*/
	bool LSIDE = (SIDE == 'L');
	int nrowA = n;
	if (LSIDE) nrowA = m;
	bool NOUNIT = (DIAG == 'N');
	bool UPPER = (UPLO == 'U');

	/*Test the input parameters.*/
	int info = 0;
	if ((!LSIDE) & (SIDE != 'R'))info = 1;
	else if ((!UPPER) & (UPLO != 'L'))info = 2;
	else if ((transA != 'N') & (transA != 'T') & (transA != 'C'))info = 3;
	else if ((DIAG != 'U') & (DIAG != 'N'))info = 4;
	else if (m < 0)info = 5;
	else if (n < 0)info = 6;
	else if (LDA < MAX(1, nrowA))info = 9;
	else if (LDB < MAX(1, m))info = 11;

	if (info != 0) {
		char name[] = "dtrmm";
		xerbla(name, info);
		return;
	}
	if (n == 0)return;
	int i, j;
	if (alpha == zero) {
		for (i = 0; i < m; i++) {
			for (j = 0; j < m; j++) {
				B[i][j] = zero;
			}
		}
		return;
	}
	int l;
	double temp;
	if (LSIDE) {
		if (transA == 'N') {
			/*Form  B := alpha*A*B.*/
			if (UPPER) {
				for (i = 0; i < m; i++) {
					if (NOUNIT)temp = alpha * A[i][i];
					else temp = alpha;
					for (j = 0; j < n; j++) {
						B[i][j] = temp * B[i][j];
					}
					for (l = i + 1; l < m; l++) {
						if (A[i][l] != 0) {
							temp = alpha * A[i][l];
							for (j = 0; j < n; j++) {
								B[i][j] = B[i][j] + B[l][j] * temp;
							}
						}
					}
				}
			}
			else {
				for (i = m - 1; i >= 0; i = i - 1) {
					temp = alpha * A[i][i];
					for (j = n - 1; j >= 0; j = j - 1) {
						B[i][j] = temp * B[i][j];
					}
					for (l = i - 1; l >= 0; l = l - 1) {
						if (A[i][l] != 0) {
							temp = alpha * A[i][l];
							for (j = n - 1; j >= 0; j = j - 1) {
								B[i][j] = B[i][j] + B[l][j] * temp;
							}
						}
					}
				}
			}
		}
		else {
			/*Form  B := alpha*A'*B*/
			if (UPPER) {
				for (i = m - 1; i >= 0; i = i - 1) {
					for (l = m - 1; l > i; l = l - 1) {
						if (A[i][l] != 0) {
							temp = alpha * A[i][l];
							for (j = n - 1; j >= 0; j = j - 1) {
								B[l][j] = B[l][j] + temp * B[i][j];
							}
						}

					}
					if (NOUNIT)temp = A[i][i] * alpha;
					else temp = alpha;
					for (j = n - 1; j >= 0; j = j - 1) {
						B[i][j] = temp * B[i][j];
					}

				}
			}
			else {
				for (i = 0; i < m; i++) {
					for (l = 0; l < i; l++) {
						if (A[i][l] != 0) {
							temp = alpha * A[i][l];
							for (j = 0; j < n; j++) {
								B[l][j] = B[l][j] + temp * B[i][j];
							}
						}
					}
					if (NOUNIT)temp = A[i][i] * alpha;
					else temp = alpha;
					for (j = 0; j < n; j++) {
						B[i][j] = temp * B[i][j];
					}
				}
			}
		}

	}
	else {
		if (transA == 'N') {
			/*Form  B := alpha*B*A.*/
			if (UPPER) {
				for (i = m - 1; i >= 0; i = i - 1) {
					for (j = n - 1; j >= 0; j = j - 1) {
						if (B[i][j] != 0) {
							temp = B[i][j] * alpha;
							for (l = n - 1; l > j; l = l - 1) {
								B[i][l] = B[i][l] + temp * A[j][l];
							}
							if (NOUNIT)  B[i][j] = temp * A[j][j];
							else B[i][j] = temp;
						}
					}
				}
			}
			else {
				for (i = 0; i < m; i++) {
					for (j = 0; j < n; j++) {
						if (B[i][j] != 0) {
							temp = B[i][j] * alpha;
							for (l = 0; l < j; l++) {
								B[i][l] = B[i][l] + temp * A[j][l];
							}
							if (NOUNIT)  B[i][j] = temp * A[j][j];
							else B[i][j] = temp;
						}
					}
				}
			}
		}
		else {
			/*Form  B := alpha*B*A'.*/
			if (UPPER) {
				for (i = 0; i < m; i++) {
					for (j = 0; j < n; j++) {
						if (NOUNIT)temp = B[i][j] * A[j][j];
						else temp = B[i][j];
						for (l = j + 1; l < n; l++) {
							temp = temp + B[i][l] * A[j][l];
						}
						B[i][j] = alpha * temp;
					}
				}
			}
			else {
				for (i = m - 1; i >= 0; i = i - 1) {
					for (j = n - 1; j >= 0; j = j - 1) {
						if (NOUNIT)temp = B[i][j] * A[j][j];
						else temp = B[i][j];
						for (l = j - 1; l >= 0; l = l - 1) {
							temp = temp + B[i][l] * A[j][l];
						}
						B[i][j] = alpha * temp;
					}
				}
			}
		}
	}
}

/*求解方程op(A)X=alpha*B,Xop(A)=alpha*B*/
void dtrsm(char SIDE, char UPLO, char transA, char DIAG, int  m, int n, double alpha,
	vector<vector<double>>& A, int  LDA, vector<vector<double>>& B, int LDB) {
	bool LSIDE = (SIDE == 'L');
	int nrowA = n;
	if (LSIDE)nrowA = m;
	bool NOUNIT = (DIAG == 'N'), UPPER = (UPLO == 'U');
	int info = 0;
	if ((!LSIDE) & (SIDE != 'R'))info = 1;
	else if ((!UPPER) & (UPLO != 'L'))info = 2;
	else if ((transA != 'N') & (transA != 'T') & (transA != 'C'))info = 3;
	else if ((DIAG != 'U') & (DIAG != 'N'))info = 4;
	else if (m < 0)info = 5;
	else if (n < 0)info = 6;
	else if (LDA < MAX(1, nrowA))info = 9;
	else if (LDB < MAX(1, m))info = 11;

	if (info != 0) {
		char name[] = "dtrsm";
		xerbla(name, info);
		return;
	}
	if (n == 0)return;
	int i, j;
	if (alpha == zero) {
		for (i = 0; i < m; i++) {
			for (j = 0; j < m; j++) {
				B[i][j] = zero;
			}
		}
		return;
	}
	int l;
	double temp;
	if (LSIDE) {
		if (transA == 'N') {
			/*Form  B := alpha*inv( A )*B.*/
			if (UPPER) {
				for (i = m - 1; i >= 0; i = i - 1) {
					if (alpha != one) {
						for (j = n - 1; j >= 0; j = j - 1) {
							B[i][j] = alpha * B[i][j];
						}
					}
					for (l = m - 1; l > i; l = l - 1) {
						if (A[i][l] != 0) {
							for (j = n - 1; j >= 0; j = j - 1) {
								B[i][j] = B[i][j] - A[i][l] * B[l][j];
							}
						}
					}
					if (NOUNIT) {
						temp = 1 / A[i][i];
						for (j = n - 1; j >= 0; j = j - 1) {
							B[i][j] = B[i][j] * temp;
						}
					}

				}
			}
			else {
				for (i = 0; i < m; i++) {
					if (alpha != one) {
						for (j = 0; j < n; j++) {
							B[i][j] = alpha * B[i][j];
						}

					}
					for (l = 0; l < i; l++) {
						if (A[i][l] != 0) {
							for (j = 0; j < n; j++) {
								B[i][j] = B[i][j] - A[i][l] * B[l][j];
							}
						}
					}
					if (NOUNIT) {
						temp = 1 / A[i][i];
						for (j = 0; j < n; j++) {
							B[i][j] = B[i][j] * temp;
						}
					}

				}

			}
		}
		else {
			/*Form  B := alpha*inv( A' )*B.*/
			if (UPPER) {
				for (i = 0; i < m; i++) {
					if (NOUNIT) {
						temp = 1 / A[i][i];
						for (j = 0; j < n; j++) {
							B[i][j] = temp * B[i][j];
						}
					}
					for (l = i + 1; l < m; l++) {
						if (A[i][l] != 0) {
							temp = A[i][l];
							for (j = 0; j < n; j++) {
								B[l][j] = B[l][j] - temp * B[i][j];
							}

						}

					}
					if (alpha != one) {
						for (j = 0; j < n; j++) {
							B[i][j] = alpha * B[i][j];
						}
					}
				}
			}
			else {
				for (i = n - 1; i >= 0; i = i - 1) {
					if (NOUNIT) {
						temp = 1 / A[i][i];
						for (j = n - 1; j >= 0; j = j - 1) {
							B[i][j] = temp * B[i][j];
						}
					}
					for (l = i - 1; l >= 0; l = l - 1) {
						if (A[i][l] != 0) {
							temp = A[i][l];
							for (j = n - 1; j >= 0; j = j - 1) {
								B[l][j] = B[l][j] - temp * B[i][j];
							}

						}

					}
					if (alpha != one) {
						for (j = n - 1; j >= 0; j = j - 1) {
							B[i][j] = alpha * B[i][j];
						}
					}
				}
			}
		}
	}
	else {
		if (transA == 'N') {
			/*Form  B := alpha*B*inv( A ).*/
			if (UPPER) {
				for (i = 0; i < m; i++) {
					if (alpha != one) {
						for (j = 0; j < n; j++) {
							B[i][j] = B[i][j] * alpha;
						}
					}
					for (j = 0; j < n; j++) {
						if (NOUNIT) {
							B[i][j] = B[i][j] / A[j][j];
						}
						for (l = j + 1; l < n; l++) {
							B[i][l] = B[i][l] - A[j][l] * B[i][j];
						}
					}
				}

			}
			else {
				for (i = m - 1; i >= 0; i = i - 1) {
					if (alpha != one) {
						for (j = n - 1; j >= 0; j = j - 1) {
							B[i][j] = B[i][j] * alpha;
						}
					}
					for (j = n - 1; j >= 0; j = j - 1) {
						if (NOUNIT) {
							B[i][j] = B[i][j] / A[j][j];
						}
						for (l = j - 1; l >= 0; l = l - 1) {
							B[i][l] = B[i][l] - A[j][l] * B[i][j];
						}
					}
				}
			}
		}
		else {
			/*Form  B := alpha*B*inv( A' ).*/
			if (UPPER) {
				for (i = m - 1; i >= 0; i = i - 1) {
					if (alpha != one) {
						for (j = n - 1; j >= 0; j = j - 1) {
							B[i][j] = B[i][j] * alpha;
						}
					}
					for (j = n - 1; j >= 0; j = j - 1) {
						for (l = n - 1; l > j; l = l - 1) {
							B[i][j] = B[i][j] - A[j][l] * B[i][l];
						}
						if (NOUNIT) {
							B[i][j] = B[i][j] / A[j][j];
						}
					}
				}
			}
			else {
				for (i = 0; i < m; i = i++) {
					if (alpha != one) {
						for (j = 0; j < n; j = j++) {
							B[i][j] = B[i][j] * alpha;
						}
					}
					for (j = 0; j < n; j++) {
						for (l = 0; l < j; l++) {
							B[i][j] = B[i][j] - A[j][l] * B[i][l];
						}
						if (NOUNIT) {
							B[i][j] = B[i][j] / A[j][j];
						}
					}
				}

			}
		}

	}

}



