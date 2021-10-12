#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "GradientLaplacian.h"

int main(int argc, char *argv[]) {
	int N = 256, M = 256, Q = 255;
	
	// test lenna image
	ImageType lenna(N, M, Q);
	readImage("../images/lenna.pgm", lenna);
	getGradient("lenna", lenna, 1);						// perform prewitt
	getGradient("lenna", lenna, 2);						// perform sobel
	getLaplacian("lenna", lenna);						// perform laplacian
	
	// test sf image
	ImageType sf(N, M, Q);
	readImage("../images/sf.pgm", sf);
	getGradient("sf", sf, 1);
	getGradient("sf", sf, 2);
	getLaplacian("sf", sf);

	return 0;
}
