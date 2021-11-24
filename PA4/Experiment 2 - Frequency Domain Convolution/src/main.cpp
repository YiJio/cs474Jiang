#include <iostream>
#include <fstream>
#include <cstdlib>

#include "image.h"
#include "fft.h"

int main(int argc, char *argv[]) {
	int N = 256, M = 256, Q = 255;
	std::complex<float>* transform = new std::complex<float>[N * M];
	std::complex<float>* before = new std::complex<float>[N * M];
	std::complex<float>* after = new std::complex<float>[N * M];
	
	// read original lenna image and get spectrum
	ImageType lenna(N, M, Q);
	readImage("../images/lenna.pgm", lenna);
	transformImage(lenna, before, 1);
	getImageSpectrum("lenna_before_spectrum", before, N, M, true);

	// perform filtering in spatial and frequency domain
	spatialFilter("lenna", lenna);
	frequencyFilter("lenna", lenna, transform);

	// read new filtered image from frequency domain and get spectrum
	ImageType lenna2(N, M, Q);
	readImage("../images/lenna_frequency.pgm", lenna2);
	transformImage(lenna2, after, 1);
	getImageSpectrum("lenna_after_spectrum", after, N, M, true);

	return 0;
}
