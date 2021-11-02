#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "fft.h"

int main(int argc, char *argv[]) {
	int N = 512, M = 512, Q = 255;
	
	// generate image 32 to test
	ImageType img32(N, M, Q);
	generateImage(32, 32, Q, 0);
	readImage("../images/img32.pgm", img32);
	std::complex<float>* transform32a = new std::complex<float>[N * M];
	shiftImage("img32", img32, transform32a, 32, 0);
	std::complex<float>* transform32b = new std::complex<float>[N * M];
	shiftImage("img32", img32, transform32b, 32, 1);
	
	// generate image 64 to test
	ImageType img64(N, M, Q);
	generateImage(64, 64, Q, 0);
	readImage("../images/img64.pgm", img64);
	std::complex<float>* transform64 = new std::complex<float>[N * M];
	shiftImage("img64", img64, transform64, 64, 0);
	shiftImage("img64", img64, transform64, 64, 1);
	
	// generate image 128 to test
	ImageType img128(N, M, Q);
	generateImage(128, 128, Q, 0);
	readImage("../images/img128.pgm", img128);
	std::complex<float>* transform128 = new std::complex<float>[N * M];
	shiftImage("img128", img128, transform128, 128, 0);
	shiftImage("img128", img128, transform128, 128, 1);

	return 0;
}
