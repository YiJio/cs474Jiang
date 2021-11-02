#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "fft.h"

int main(int argc, char *argv[]) {
	int N = 256, M = 256, Q = 255;
	
	// generate lenna to test
	ImageType lenna(N, M, Q);
	readImage("../images/lenna.pgm", lenna);
	std::complex<float>* transform = new std::complex<float>[N * M];
	shiftImage("lenna", lenna, transform, 0);
	
	/*
	ImageType lennaNew(N, M, Q);
	
	std::complex<float> test[N * M];
	std::copy(transform, transform + (N * M), test);
	for(int i = 0; i < N * M; i++) {
		//test[i] = {std::abs(transform[i]), 0};
		double theta = atan2(transform[i].imag(), transform[i].real());
		test[i] = {cos(theta), sin(theta)};
		//transform[i].std::abs(transform[i]);
		//transform[i] /= std::abs(transform[i]);
	}
	
	fft2D(test, N, M, 1); int value;
	float min = std::round(test[0].real()), max = std::round(test[0].real());	

	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			//float real = std::round(test[i * M + j].real());
			//lennaNew.setPixelVal(i, j, real);
			float real = std::abs(test[i * M + j]);
			lennaNew.setPixelVal(i, j, real);
			min = std::min(min, real);
			max = std::max(max, real);
		}
	}
	
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			lennaNew.getPixelVal(i, j, value);
			int newvalue = 255l * (value - min) / (max - min);
			lennaNew.setPixelVal(i, j, newvalue);
		}
	}
	
	std::string fname = "../images/lenna_mag1.pgm";
	char *imageFile = new char[fname.length() + 1];
	strcpy(imageFile, fname.c_str());
	writeImage(imageFile, lennaNew);
	delete[] imageFile;*/
	
	return 0;
}
