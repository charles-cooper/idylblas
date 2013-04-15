#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>

int
main()
{
	unsigned int const len = 1000;
	double x[len];
	double y[len];
	for (int i = 0; i < len; i++)
	{
		x[i] = rand() / (double)RAND_MAX;
		y[i] = rand() / (double)RAND_MAX;
	}
	double ret;
	int inc = 1;
	for (int i = 0; i < len*len; i++)
		ret = cblas_ddot(len, x, inc, y, inc);
	return (int)ret;
}

