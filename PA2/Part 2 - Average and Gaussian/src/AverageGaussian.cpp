#include "AverageGaussian.h"

/**
 * This function computes for average smooth filtering of the image. Based on the mask size being
 * passed in, it will create a window that is of that mask size, which helps to determine the offset
 * and bounds of the window. It will make a window and then take the accumulated sum from all the
 * window's values. The sum is then averaged by the number of pixels in mask and returned.
 * @param: image ImageType reference, integer values for original row and col for current pixel,
 * integer for max bound
 * @return: integer average
 */
int computeAverage(ImageType& image, int row, int col, int masksize, int max) {
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
	
	// get sums of average and average it
	int sum = 0;
	for(int i = 0; i < masksize; i++) {
		for(int j = 0; j < masksize; j++) {
			sum += window[i][j];
		}
	}
	int avg = sum / (masksize * masksize);
	
	return avg;
}

/**
 * This function takes in the image's info and then calls computeAverage to compute the average of the
 * current pixel's window and writes the new image out.
 * @param: fname character array to write image, image ImageType reference, integer mask size
 * @pre: all original values in variables
 * @post: average image written
 * @return: none
 */
void getAverage(char fname[], ImageType& image, int masksize) {
	// variables
	int M, N, Q, value;
	image.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);

	// perform average
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			value = computeAverage(image, i, j, masksize, N);
			newImage.setPixelVal(i, j, value);
		}
	}

	// average file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_average_mask" + std::to_string(masksize) + ".pgm";
	char *newImageFile = new char[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
	delete[] newImageFile;
}

/**
 * This function computes for gaussian smooth filtering of the image. Based on the mask size being
 * passed in, it will create a window that is of that mask size, which helps to determine the offset
 * and bounds of the window. It will make a window and then take the accumulated sum from all the
 * window's values multiplied by the gaussian mask value. This value is determined by the mask size,
 * which uses either mask7 or mask15. The sum is then averaged by the sum of all weights in the mask
 * and returned.
 * @param: image ImageType reference, integer values for original row and col for current pixel,
 * integer for max bound
 * @return: integer average
 */
int computeGaussian(ImageType& image, int row, int col, int masksize, int max) {	
	// gaussian mask 7x7
	int mask7[7][7] = {
		{1,1,2,2,2,1,1},
		{1,2,2,4,2,2,1},
		{2,2,4,8,4,2,2},
		{2,4,8,16,8,4,2},
		{2,2,4,8,4,2,2},
		{1,2,2,4,2,2,1},
		{1,1,2,2,2,1,1}
	};
	// gaussian mask 15x15
	int mask15[15][15] = {
		{2,2,3,4,5,5,6,6,6,5,5,4,3,2,2},
		{2,3,4,5,7,7,8,8,8,7,7,5,4,3,2},
		{3,4,6,7,9,10,10,11,10,10,9,7,6,4,3},
		{4,5,7,9,10,12,13,13,13,12,10,9,7,5,4},
		{5,7,9,11,13,14,15,16,15,14,13,11,9,7,5},
		{5,7,10,12,14,16,17,18,17,16,14,12,10,7, 5},
		{6,8,10,13,15,17,19,19,19,17,15,13,10,8,6},
		{6,8,11,13,16,18,19,20,19,18,16,13,11,8,6},
		{6,8,10,13,15,17,19,19,19,17,15,13,10,8,6},
		{5,7,10,12,14,16,17,18,17,16,14,12,10,7,5},
		{5,7,9,11,13,14,15,16,15,14,13,11,9,7,5},
		{4,5,7,9,10,12,13,13,13,12,10,9,7,5,4},
		{3,4,6,7,9,10,10,11,10,10,9,7,6,4,3},
		{2,3,4,5,7,7,8,8,8,7,7,5,4,3,2},
		{2,2,3,4,5,5,6,6,6,5,5,4,3,2,2}
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
			if(masksize == 7) {
				sum += window[i][j] * mask7[i][j];
				summask += mask7[i][j];
			}
			else if(masksize == 15) {
				sum += window[i][j] * mask15[i][j];
				summask += mask15[i][j];
			}			
		}
	}
	int avg = sum / summask;
	
	return avg;
}

/**
 * This function takes in the image's info and then calls computeGaussian to compute the gaussian filter
 * of the current pixel's window and writes the new image out.
 * @param: fname character array to write image, image ImageType reference, integer mask size
 * @pre: all original values in variables
 * @post: gaussian image written
 * @return: none
 */
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
	char *newImageFile = new char[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
	delete[] newImageFile;
}
