#include "UnsharpHighBoost.h"

int computeGaussian(ImageType& image, int row, int col, int masksize, int max) {
	// gaussian mask 7x7
	int mask[7][7] = {
		{1,1,2,2,2,1,1},
		{1,2,2,4,2,2,1},
		{2,2,4,8,4,2,2},
		{2,4,8,16,8,4,2},
		{2,2,4,8,4,2,2},
		{1,2,2,4,2,2,1},
		{1,1,2,2,2,1,1}
	};

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
	
	// grab sums of gaussian and return average of it by summask
	int sum = 0;
	int summask = 0;
	for(int i = 0; i < masksize; i++) {
		for(int j = 0; j < masksize; j++) {
			sum += window[i][j] * mask[i][j];
			summask += mask[i][j];
		}
	}
	int avg = sum / summask;
	
	return avg;
}

void getGaussian(char fname[], ImageType& image, int masksize) {
	// variables
	int M, N, Q, value;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);

	// perform gaussian
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			value = computeGaussian(image, i, j, masksize, N);
			newImage.setPixelVal(i, j, value);
		}
	}

	// gaussian file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_gaussian_mask" + std::to_string(masksize) + ".pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
}

void getSharp(char fname[], ImageType& image, ImageType& mask, int k) {
	// variables
	int M, N, Q, value, valueM;
	image.getImageInfo(N, M, Q);
	mask.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);
	
	// add image with weighted mask
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			image.getPixelVal(i, j, value);
			mask.getPixelVal(i, j, valueM);
			newImage.setPixelVal(i, j, value + (k*valueM));
		}
	}
	
	// sharp filtering file names and writing
	std::string v;
	if(k == 1) { v = "="; }
	else { v = ">"; }
	std::string newfname = "../images/" + std::string(fname) + "_sharp_k" + v + "1.pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
}

void getSharpMask(char fname[], ImageType& image, ImageType& lpImage) {
	// variables
	int M, N, Q, value, valueLP;
	image.getImageInfo(N, M, Q);
	lpImage.getImageInfo(N, M, Q);
	ImageType maskImage(N, M, Q);
	
	// subtract image from low-pass image to get mask
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			image.getPixelVal(i, j, value);
			lpImage.getPixelVal(i, j, valueLP);
			maskImage.setPixelVal(i, j, value - valueLP);
		}
	}
	
	// sharp mask file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_sharp_mask.pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, maskImage);
}
