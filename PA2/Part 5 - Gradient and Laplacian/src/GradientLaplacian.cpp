#include "GradientLaplacian.h"

// prewitt fy mask
int fyPrewitt[3][3] = {
	{-1,-1,-1},
	{0,0,0},
	{1,1,1}
};
// prewitt fx mask
int fxPrewitt[3][3] = {
	{-1,0,1},
	{-1,0,1},
	{-1,0,1}
};
// sobel fy mask
int fySobel[3][3] = {
	{-1,-2,-1},
	{0,0,0},
	{1,2,1}
};
// sobel fx mask
int fxSobel[3][3] = {
	{-1,0,1},
	{-2,0,2},
	{-1,0,1}
};
// laplacian mask
int laplacian[3][3] = {
	{0,1,0},
	{1,-4,1},
	{0,1,0}
};

/**
 * This function computes for gradient sharp filtering of the image. Based on the mask being passed in,
 * it will take either the prewitt or sobel partial derivative x or y masks, which uses a mask size of
 * 3, making the window be 3x3. The function also takes partial derivatives x, y, and the gradient
 * magnitude references. These accumulated sums are computed by multiplying the window and the desired
 * mask value, or computing the square root of both x and y squared for magnitude.
 * @param: image ImageType reference, integer values for original row and col for current pixel,
 * integer for max bound, integer mask size, integer x, y, mag references
 * @pre: original x, y, mag values from references
 * @post: new computed x, y, mag values to references
 * @return: none
 */
void computeGradient(ImageType& image, int row, int col, int max, int mask, int& x, int& y, int& mag) {	
	int masksize = 3;
	
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
	
	// compute partial derivatives for x, y depending on mask
	int fxi = 0, fyi = 0;
	for(int i = 0; i < masksize; i++) {
		for(int j = 0; j < masksize; j++) {
			if(mask == 1) {
				fxi += window[i][j] * fxPrewitt[i][j];
				fyi += window[i][j] * fyPrewitt[i][j];
			}
			else if(mask == 2) {
				fxi += window[i][j] * fxSobel[i][j];
				fyi += window[i][j] * fySobel[i][j];			
			}
		}
	}
	
	// calculate gradient magnitude
	mag = std::sqrt(std::pow(fxi, 2) + std::pow(fyi, 2));
	x = fxi;
	y = fyi;
}

/**
 * This function takes in the image's info, attempts to generate x, y, and magnitude images, calls
 * computeGradient based on mask and mode to store certain values, takes the normalization values of
 * x and y, and writes three different images out.
 * @param: fname character array to write image, image ImageType reference, integer mask preference
 * @pre: all original values in variables
 * @post: partial derivatives x, y, and gradient magnitude images written
 * @return: none
 */
void getGradient(char fname[], ImageType& image, int mask) {
	// variables
	int M, N, Q, valueX = 0, valueY = 0, valueMag = 0;
	image.getImageInfo(N, M, Q);
	ImageType xImage(N, M, Q);
	ImageType yImage(N, M, Q);
	ImageType magImage(N, M, Q);
	
	// normalization map min max
	int minX, maxX, minY, maxY;
	image.getPixelVal(0, 0, minX);
	minY = maxY = maxX = minX;

	// perform gradient
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			// compute values for x, y, magnitude by passing mask
			computeGradient(image, i, j, N, mask, valueX, valueY, valueMag);
			if(valueX < minX) { minX = valueX; }
			if(valueX > maxX) { maxX = valueX; }
			if(valueY < minY) { minY = valueY; }
			if(valueY > maxY) { maxY = valueY; }
			// set values to images based on x, y, magnitude
			xImage.setPixelVal(i, j, valueX);
			yImage.setPixelVal(i, j, valueY);
			magImage.setPixelVal(i, j, valueMag);
		}
	}
	
	// normalize to range 0-255
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			// for partial derivative x
			xImage.getPixelVal(i, j, valueX);
			int newvalueX = 255l * (valueX - minX) / (maxX - minX);
			xImage.setPixelVal(i, j, newvalueX);
			// for partial derivative y
			yImage.getPixelVal(i, j, valueY);
			int newvalueY = 255l * (valueY - minY) / (maxY - minY);
			yImage.setPixelVal(i, j, newvalueY);
		}
	}

	// gradient file names and writing
	std::string v;
	if(mask == 1) { v = "prewitt"; }
	else if(mask == 2) { v = "sobel"; }
	
	// for partial derivative x
	std::string newxname = "../images/" + std::string(fname) + "_" + v + "X.pgm";
	char *newXFile = new char[newxname.length() + 1];
	strcpy(newXFile, newxname.c_str());
	writeImage(newXFile, xImage);
	delete[] newXFile;
	
	// for partial derivative y
	std::string newyname = "../images/" + std::string(fname) + "_" + v + "Y.pgm";
	char *newYFile = new char[newyname.length() + 1];
	strcpy(newYFile, newyname.c_str());
	writeImage(newYFile, yImage);
	delete[] newYFile;
	
	// for gradient magnitude
	std::string newmagname = "../images/" + std::string(fname) + "_" + v + "Mag.pgm";
	char *newMagFile = new char[newyname.length() + 1];
	strcpy(newMagFile, newmagname.c_str());
	writeImage(newMagFile, magImage);
	delete[] newMagFile;
}

/**
 * This function computes for laplacian sharp filtering of the image. The function uses the laplacian
 * mask of 3x3, making the window be 3x3. It computes for the sum of all window values multiplied by
 * the desired mask value, which is then returned for the current pixel.
 * @param: image ImageType reference, integer values for original row and col for current pixel,
 * integer for max bound
 * @return: integer sum
 */
int computeLaplacian(ImageType& image, int row, int col, int max) {
	int masksize = 3;
	
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
	
	// get sums of window and laplacian mask
	int sum = 0;
	for(int i = 0; i < masksize; i++) {
		for(int j = 0; j < masksize; j++) {
			sum += window[i][j] * laplacian[i][j];
		}
	}
	
	return sum;
}

/**
 * This function takes in the image's info and then calls computeLaplacian to compute the laplacian filter
 * of the current pixel's window and writes the new image out.
 * @param: fname character array to write image, image ImageType reference
 * @pre: all original values in variables
 * @post: laplacian image written
 * @return: none
 */
void getLaplacian(char fname[], ImageType& image) {
	// variables
	int M, N, Q, value;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);

	// normalization map min max	
	int min, max;
	image.getPixelVal(0, 0, min);
	max = min;

	// perform laplacian
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			value = computeLaplacian(image, i, j, N);
			if(value < min) { min = value; }
			if(value > max) { max = value; }
			newImage.setPixelVal(i, j, value);
		}
	}
	
	// normalize to range 0-255
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			newImage.getPixelVal(i, j, value);
			int newvalue = 255l * (value - min) / (max - min);
			newImage.setPixelVal(i, j, newvalue);
		}
	}

	// laplacian file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_laplacian.pgm";
	char *newImageFile = new char[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
	delete[] newImageFile;
}
