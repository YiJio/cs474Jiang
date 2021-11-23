#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "NoiseRemoval.h"

int main(int argc, char *argv[]) {
	int N = 512, M = 512, Q = 255;
	
	// test with boy image
	ImageType boy(N, M, Q);
	readImage("../images/boy_noisy.pgm", boy);
	std::complex<float>* transform = new std::complex<float>[N * M];
	// use band-reject filter: boy image, boy image spectrum, filtered image, filtered spectrum
	transformImage("boy_noisy", boy, transform);
	//removeNoise(boy, transform, 0);
	// use notch filter: boy image, boy image spectrum, filtered image, filtered spectrum
	//removeNoise(boy, transform, 1);
	// use gaussian filtering 15x15
	getGaussian("boy_noisy", boy, 7);
	getGaussian("boy_noisy", boy, 15);
	

	return 0;
}
