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

	// get before image spectrum
	transformImage(boy, transform, 1);
	getImage("boy_noisy_spectrum", transform, N, M, true, 1);
	
	// perform band-reject filter (using Butterworth by default)
	computeBand("boy_noisy", boy, 1, 2, 35);
	
	// get after image spectrum using band-reject filter
	ImageType boyband(N, M, Q);
	readImage("../images/boy_noisy_band.pgm", boyband);
	std::complex<float>* band = new std::complex<float>[N * M];
	transformImage(boyband, band, 1);
	getImage("boy_noisy_spectrum_band", band, N, M, true, 1);
	
	// perform notch-reject filter (using Butterworth by default)
	computeNotch("boy_noisy", boy, 1, 5, 2, 16, 32);
	
	// get after image spectrum using notch-reject filter
	ImageType boynotch(N, M, Q);
	readImage("../images/boy_noisy_notch.pgm", boynotch);
	std::complex<float>* notch = new std::complex<float>[N * M];
	transformImage(boynotch, notch, 1);
	getImage("boy_noisy_spectrum_notch", notch, N, M, true, 1);
	
	// use gaussian filtering
	getGaussian("boy_noisy", boy, 7);
	getGaussian("boy_noisy", boy, 15);	

	return 0;
}
