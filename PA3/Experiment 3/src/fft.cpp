#include "fft.h"
#include <omp.h>

void fft(std::complex<float> data[], int n, int isign, int r) {
	// variables
	int j = n / 2;
	
	for(int i = 1; i < (n - 1); i++) {
		if(j > i) { std::swap(data[i * r], data[j * r]); }
		int m = n / 2;
		while(m > 0 && j & m) {
			j ^= m;
			m >>= 1;
		}
		j |= m;
	}
	
	for(int M = 2; M <= n; M *= 2) {
		double theta = isign * (2 * M_PI / M);
		std::complex<double> W_M = {cos(theta), sin(theta)};
		std::complex<double> W_Mu = {1, 0};
		for(int u = 0; u < M / 2; u++) {
			for(int i = u; i < n; i += M) {
				int j = i + M / 2;
				std::complex<double> curr = std::complex<double>(data[j * r]) * W_Mu;
				data[j * r] = std::complex<double>(data[i * r]) - curr;
				data[i * r] += curr;
			}
			W_Mu *= W_M;
		}
	}
}

void fft2D(std::complex<float> data[], int N, int M, int isign) {
	float factor = 1 / sqrt(N * M);
	#pragma omp parallel
	{
		#pragma omp for
		for(int i = 0; i < N; i++) {
			fft(data + i * M, M, isign);
		}
		#pragma omp for
		for(int i = 0; i < M; i++) {
			fft(data + i, N, isign, M);
		}
		#pragma omp for
		for(int i = 0; i < N * M; i++) {
			data[i] *= factor;
		}
	}
	
}

void generateImage(int n, int m, int Q, int mode) {
	int N = 512, M = 512;
	ImageType image(N, M, Q);
	
	// fill image black with nxn center white
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			if((i>(N/2-n/2)) && (i<=(N-N/2+n/2)) && (j>(M/2-m/2)) && (j<=(M-M/2+m/2))) {
				image.setPixelVal(i, j, 255);
			} else { image.setPixelVal(i, j, 0); }
		}
	}

	// file names and writing
	std::string v = "";
	if(mode == 1) { v = "_mag_shift_n"; }
	else if(mode == 2) { v = "_mag_shift_y"; }
	else if(mode == 3) { v = "mag_shift_yLog"; }
	std::string fname = "../images/img" + std::to_string(n) + v + ".pgm";
	char *imageFile = new char[fname.length() + 1];
	strcpy(imageFile, fname.c_str());
	writeImage(imageFile, image);
	delete[] imageFile;
}

void shiftImage(char fname[], ImageType& image, std::complex<float> transform[], int mode) {
	// variables
	int M, N, Q, value, curr;
	image.getImageInfo(N, M, Q);
	
	// perform shift if mode is 1
	if(mode == 1) {
		#pragma omp parallel for collapse(2)
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < M; j++) {
				// translate/shift magnitude to center of frequency domain
				if((i + j) % 2 == 0) { curr = 1; }
				else { curr = -1; }
				image.getPixelVal(i, j, value);
				transform[i * M + j] = {value * curr, 0};
			}
		}
	} else {
		// keep current values
		#pragma omp parallel for
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < M; j++) {
				image.getPixelVal(i, j, value);
				transform[i * M + j] = {value, 0};
			}
		}
	}
	fft2D(transform, N, M, -1);
	computeImage(transform, N, M, 0);
	computeImage(transform, N, M, 1);
	computeImage(transform, N, M, 2);
}

void computeImage(std::complex<float> transform[], int N, int M, int mode) {
	// variables
	int Q = 255, value;
	ImageType image(N, M, Q);
	std::complex<float> test[N * M];
	std::copy(transform, transform + (N * M), test);

	// phase 0, mag 1, mag log 2
	for(int i = 0; i < N * M; i++) {
		//if(mode == 1 || mode == 2) { test[i] = std::abs(test[i]); }
		//else { test[i] /= std::abs(test[i]); }
		if(mode == 0 || mode == 2) { test[i] = {std::abs(transform[i]), 0}; }
		else {
			double theta = atan2(transform[i].imag(), transform[i].real());
			test[i] = {cos(theta), sin(theta)};
		}
	}
	
	fft2D(test, N, M, 1);
	float min = std::round(test[0].real()), max = std::round(test[0].real());	

	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			float real = std::round(test[i * M + j].real());
			if(mode == 2) { real = 2000 * log(1 + real); }
			image.setPixelVal(i, j, real);
			min = std::min(min, real);
			max = std::max(max, real);
		}
	}
	
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			image.getPixelVal(i, j, value);
			int newvalue = 255l * (value - min) / (max - min);
			image.setPixelVal(i, j, newvalue);
		}
	}
	
	std::string v = "";
	if(mode == 0) { v = "removePhase"; }
	else if(mode == 1) { v = "removeMag"; }
	else if(mode == 2) { v = "removeMagLog"; }
	std::string fname = "../images/lenna_" + v + ".pgm";
	char *imageFile = new char[fname.length() + 1];
	strcpy(imageFile, fname.c_str());
	writeImage(imageFile, image);
	delete[] imageFile;
}
