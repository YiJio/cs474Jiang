#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "fft.h"

int main(int argc, char *argv[]) {
	int N = 512, M = 512, Q = 255;
	
	// generate lenna image to test
	ImageType girl(N, M, Q);
	std::complex<float>* t = new std::complex<float>[N * M];
	readImage("../images/girl.pgm", girl);
	transformImage(girl, t, 1);
	getImage("girl_spectrum", t, N, M, true, 1);
	
	// apply homomorphic filtering on different gamma values
	homomorphicFilter("girl", girl, 1.5, 0.5);
	homomorphicFilter("girl", girl, 1.0, 0.0);
	homomorphicFilter("girl", girl, 1.0, 1.0);
	homomorphicFilter("girl", girl, 2.0, 0.0);
	homomorphicFilter("girl", girl, 2.0, 1.0);
	
	ImageType girl1(N, M, Q);
	std::complex<float>* t1 = new std::complex<float>[N * M];
	
	readImage("../images/girl_H15_L05.pgm", girl1);
	transformImage(girl1, t1, 1);
	getImage("girl_H15_L05_spectrum", t1, N, M, true, 1);
	
	readImage("../images/girl_H1_L0.pgm", girl1);
	transformImage(girl1, t1, 1);
	getImage("girl_H1_L0_spectrum", t1, N, M, true, 1);
	
	readImage("../images/girl_H2_L0.pgm", girl1);
	transformImage(girl1, t1, 1);
	getImage("girl_H2_L0_spectrum", t1, N, M, true, 1);
	
	return 0;
}
