#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "fft.h"

int main(int argc, char *argv[]) {
	int N = 256, M = 256, Q = 255;
	
	// generate lenna image to test
	ImageType lenna(N, M, Q);
	readImage("../images/lenna.pgm", lenna);
	std::complex<float>* transform = new std::complex<float>[N * M];
	transformImage("lenna", lenna, transform);
	
	return 0;
}
