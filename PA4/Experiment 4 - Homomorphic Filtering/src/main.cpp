#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "fft.h"

int main(int argc, char *argv[]) {
	int N = 512, M = 512, Q = 255;
	
	// generate lenna image to test
	ImageType girl(N, M, Q);
	readImage("../images/girl.pgm", girl);
	
	// apply homomorphic filtering on different gamma values
	homomorphicFilter("girl", girl, 1.5, 0.5);
	homomorphicFilter("girl", girl, 1.0, 0.0);
	homomorphicFilter("girl", girl, 2.0, 1.0);
	homomorphicFilter("girl", girl, 1.0, 1.0);
	
	return 0;
}
