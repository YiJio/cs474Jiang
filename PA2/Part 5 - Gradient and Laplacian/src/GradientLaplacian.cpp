#include "GradientLaplacian.h"

int computePrewitt(ImageType& image, int row, int col, int max) {	
	int masksize = 3;
	int fy[masksize][masksize] = {
		{-1,-1,-1},
		{0,0,0},
		{1,1,1}
	};
	int fx[masksize][masksize] = {
		{-1,0,1},
		{-1,0,1},
		{-1,0,1}
	};
	
	int k = row;
	int l = col;
	int o = masksize / 2;
	int window[masksize][masksize];
	
	for(int i = 0; i < masksize; i++, k++) {
		l = col;
		for(int j = 0; j < masksize; j++, l++) {
			if((k-o>=0 && k-o<max) && (l-o>=0 && l-o<max)) {
				image.getPixelVal(k-o, l-o, window[i][j]);
			}
			else { window[i][j] = 0; }
		}
	}
	int fxi = 0, fyi = 0;
	for(int i = 0; i < masksize; i++) {
		for(int j = 0; j < masksize; j++) {
			fxi += window[i][j] * fx[i][j];
			fyi += window[i][j] * fy[i][j];
		}
	}
	int mag = std::sqrt(std::pow(fxi, 2) + std::pow(fyi, 2));
	return mag;
}

void getPrewitt(char fname[], ImageType& image) {
	// variables
	int M, N, Q, value;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);

	// perform prewitt
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			value = computePrewitt(image, i, j, N);
			newImage.setPixelVal(i, j, value);
		}
	}

	// prewitt file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_prewitt.pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
}

int computeSobel(ImageType& image, int row, int col, int max) {
	int masksize = 3;
	int fy[masksize][masksize] = {
		{-1,-2,-1},
		{0,0,0},
		{1,2,1}
	};
	int fx[masksize][masksize] = {
		{-1,0,1},
		{-2,0,2},
		{-1,0,1}
	};	
	
	int k = row;
	int l = col;
	int o = masksize / 2;
	int window[masksize][masksize];
	
	for(int i = 0; i < masksize; i++, k++) {
		l = col;
		for(int j = 0; j < masksize; j++, l++) {
			if((k-o>=0 && k-o<max) && (l-o>=0 && l-o<max)) {
				image.getPixelVal(k-o, l-o, window[i][j]);
			}
			else { window[i][j] = 0; }
		}
	}
	int sum = 0;
	for(int i = 0; i < masksize; i++) {
		for(int j = 0; j < masksize; j++) {
			sum += window[i][j];
		}
	}
	int avg = sum / (masksize*masksize);
	return avg;
}

void getSobel(char fname[], ImageType& image) {
	// variables
	int M, N, Q, value;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);

	// perform sobel
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			value = computeSobel(image, i, j, N);
			newImage.setPixelVal(i, j, value);
		}
	}

	// sobel file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_sobel.pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
}

int computeLaplacian(ImageType& image, int row, int col, int max) {
	int masksize = 3;
	int k = row;
	int l = col;
	int o = masksize / 2;
	int window[masksize][masksize];
	
	for(int i = 0; i < masksize; i++, k++) {
		l = col;
		for(int j = 0; j < masksize; j++, l++) {
			if((k-o>=0 && k-o<max) && (l-o>=0 && l-o<max)) {
				image.getPixelVal(k-o, l-o, window[i][j]);
			}
			else { window[i][j] = 0; }
		}
	}
	
	int sum = 0;
	int summask = 0;
	for(int i = 0; i < masksize; i++) {
		for(int j = 0; j < masksize; j++) {
			
		}
	}
	int avg = sum / summask;
	
	return avg;
}

void getLaplacian(char fname[], ImageType& image) {
	// variables
	int M, N, Q, value;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);

	// perform laplacian
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			value = computeLaplacian(image, i, j, N);
			newImage.setPixelVal(i, j, value);
		}
	}

	// laplacian file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_laplacian.pgm";
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
}
