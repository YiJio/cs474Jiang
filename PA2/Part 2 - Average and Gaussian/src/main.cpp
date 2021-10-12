#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "AverageGaussian.h"

int main(int argc, char *argv[]) {
	int N = 256, M = 256, Q = 255;
	
	// test lenna image
	ImageType lenna(N, M, Q);
	readImage("../images/lenna.pgm", lenna);
	getAverage("lenna", lenna, 7);						// perform average mask 7x7
	getAverage("lenna", lenna, 15);						// perform average mask 15x15
	getGaussian("lenna", lenna, 7);						// perform gaussian mask 7x7
	getGaussian("lenna", lenna, 15);					// perform gaussian mask 15x15
	
	// test sf image
	ImageType sf(N, M, Q);
	readImage("../images/sf.pgm", sf);
	getAverage("sf", sf, 7);
	getAverage("sf", sf, 15);
	getGaussian("sf", sf, 7);
	getGaussian("sf", sf, 15);

	return 0;
}
