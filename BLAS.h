#include <math.h>
#include <vector>
#include "xerbal.h"
using namespace std;
#define zero 0.00
#define one 1.00


//求欧几里得2-范数
double dnrm2(int n, double* x, int  incx) {
	double dnrm2 = 0;
	if (n <= 0) return dnrm2;
	int i = 0, ix, j, next = 30;
	double sum = 0.00, xmax, temp, cutlo = 1e-308, cuthi = 1e308, hitest;

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
