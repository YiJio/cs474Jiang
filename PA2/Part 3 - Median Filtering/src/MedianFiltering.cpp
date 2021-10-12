#include "MedianFiltering.h"

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
 * This function computes for median smooth filtering of the image, which is of non-linear type. Based on
 * the mask size being passed in, it will create a window that is of that mask size, which helps to
 * determine the offset and bounds of the window. It will make a window, sort all values within that
 * neighborhood, and then take the median and return it.
 * @param: image ImageType reference, integer values for original row and col for current pixel,
 * integer for max bound
 * @return: integer median
 */
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
	int median = neighborhood[n / 2];
	if(n % 2 == 0) { median = (neighborhood[n / 2] + neighborhood[(n / 2) - 1]) / 2; }
	
	return median;
}

/**
 * This function takes in the image's info and then calls computeMedian to compute the median of the
 * current pixel's window and writes the new image out.
 * @param: fname character array to write image, image ImageType reference, integer mask size
 * @pre: all original values in variables
 * @post: median image written
 * @return: none
 */
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
	char *newImageFile = new char[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
	delete[] newImageFile;
}

/**
 * This function takes in the image's info and creates a copy of the image to be "salted and pepper".
 * For the full size of the original image, it will do a random chance based on variable X. It will
 * generate a random row and col location and then assign randomness of either black or white to that
 * pixel. The new corrupted image is then written out.
 * @param: fname character array to write image, image ImageType reference, integer X for
 * random chance
 * @pre: all original values in variables
 * @post: salted pepper image written
 * @return: none
 */
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
	for(int i = 0; i < N * M; ++i) {
		int chance = rand() % 100;
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
	char *newImageFile = new char[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, spImage);
	delete[] newImageFile;
}
