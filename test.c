#include  <stdio.h>
#include  <vector>
#include "myprint.h"
#include "BLAS.h"
using namespace std;

int main()
{
	double arrA[7] = {3,4,7,5,3,1e-30,1e30};
	int inca = 1;
	//incb = 1;
	double dnorm = dnrm2(7,arrA,inca);
	printArray(&dnorm,1);
	return 0;
}
