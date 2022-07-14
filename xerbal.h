#include<iostream>
using namespace std;

void xerbla(char srname[], int info);
int MAX(int x, int y);
int MIN(int x, int y);

void xerbla(char srname[], int info) {
	cout << "On entry to parameter number  had an illegal value " <<endl; 
	cout << srname << endl;
	cout << info << endl;
	return;
}

int MAX(int x, int y) {
	int max = x;
	if (y > x)max = y;
	return max;
}

int MIN(int x, int y) {
	int min = x;
	if (y < x)min = y;
	return min;
}
