#include "MedianFiltering.h"

int computeMedian(ImageType& image, int row, int col, int masksize, int max) {
	// variables
	int k = row;
	int l = col;
	int o = masksize / 2;
	int window[masksize][masksize];
	
	// get window for current pixel
	for(int i = 0; i < masksize; i++, k++) {
		l = col;
		for(int j = 0; j < masksize; j++, l++) {
			if((k-o>=0 && k-o<max) && (l-o>=0 && l-o<max)) {
				image.getPixelVal(k-o, l-o, window[i][j]);
			}
			else { window[i][j] = 0; }
		}
	}
	
	// make neighborhood to 1D array
	int x = 0;
	int neighborhood[masksize*masksize];
	for(int i = 0; i < masksize; i++) {
		for(int j = 0; j < masksize; j++) {
			neighborhood[x++] = window[i][j];
		}
	}
	
	// sort neighborhood array
	int n = sizeof(neighborhood) / sizeof(neighborhood[0]);
	std::sort(neighborhood, neighborhood + n);
	
	// take median of odd array first, then change if array is even
	int median = neighborhood[n/2];
	if(n % 2 == 0) { std::cout << "es"; median = (neighborhood[n/2] + neighborhood[(n/2)-1])/2; }
	
	return median;
}

void getMedian(char fname[], ImageType& image, int masksize) {
	// variables
	int M, N, Q, median;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);
	
	// perform median on salted pepper image
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			median = computeMedian(image, i, j, masksize, N);
			newImage.setPixelVal(i, j, median);
		}
	}
	
	// median file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_median_mask" + std::to_string(masksize) + ".pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
}

void saltedPepper(char fname[], ImageType &image, int X) {
	// variables
	int M, N, Q, value;
	image.getImageInfo(N, M, Q);
	ImageType spImage(N, M, Q);
	
	// create copy image to be corrupted
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			image.getPixelVal(i, j, value);
			spImage.setPixelVal(i, j, value);
		}
	}
	
	// perform salt and pepper, corrupt copy image
	for(int i = 0; i < N*M; ++i) {
		int chance = random();
		chance *= 100;
		// if under X%, do corruption stuff
		if(chance < X) {
			// randomly pick pixel location
			int x = rand() % 256;
			int y = rand() % 256;
			// randomly assign black or white on pixel
			int bw = rand() % 2;
			if(bw == 0) { spImage.setPixelVal(x, y, 0); }
			else if(bw == 1) { spImage.setPixelVal(x, y, 255); }
		}
	}
	
	// salted pepper file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_X" + std::to_string(X) + ".pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, spImage);
}
