#include "UnsharpHighBoost.h"

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

/**
 * This function takes in the image and mask so that it can compute the new sharpened image based on
 * value of k. It computes the function of new_image = original_image + k * mask. The sharpened image
 * is then normalized, as sharpening has negative values out of the 0-255 range. The sharpened image
 * then gets written out.
 * @param: fname character array to write image, image ImageType reference, mask ImageType reference,
 * integer k amount/weight
 * @pre: all original values in variables
 * @post: sharpened image written
 * @return: none
 */
void getSharp(char fname[], ImageType& image, ImageType& mask, int k) {
	// variables
	int M, N, Q, value, valueM;
	image.getImageInfo(N, M, Q);
	mask.getImageInfo(N, M, Q);
	ImageType newImage(N, M, Q);
	
	// normalization map min max
	int min, max;
	image.getPixelVal(0, 0, min);
	max = min;

	// add image with weighted mask
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			image.getPixelVal(i, j, value);
			mask.getPixelVal(i, j, valueM);
			int curr = value + (k * valueM);
			if(curr < min) { min = curr; }
			if(curr > max) { max = curr; }
			newImage.setPixelVal(i, j, curr);
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
	
	// sharp filtering file names and writing
	std::string v;
	if(k == 1) { v = "="; }
	else { v = ">"; }
	std::string newfname = "../images/" + std::string(fname) + "_sharp_k" + v + "1.pgm";
	char *newImageFile = new char[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
	delete[] newImageFile;
}

/**
 * This function takes in the image and low-pass (smoothed) image so that it can compute for the sharp
 * mask. It computes the function of mask = original_image - lp_image. The sharp mask is then normalized,
 * as the subtraction could lead to negative values out of the 0-255 range. The sharp mask then gets
 * written out.
 * @param: fname character array to write image, image ImageType reference, lpImage ImageType reference
 * @pre: all original values in variables
 * @post: sharp mask written
 * @return: none
 */
void getSharpMask(char fname[], ImageType& image, ImageType& lpImage) {
	// variables
	int M, N, Q, value, valueLP;
	image.getImageInfo(N, M, Q);
	lpImage.getImageInfo(N, M, Q);
	ImageType maskImage(N, M, Q);

	// normalization map min max	
	int min, max;
	image.getPixelVal(0, 0, min);
	max = min;

	// subtract image from low-pass image to get mask
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			image.getPixelVal(i, j, value);
			lpImage.getPixelVal(i, j, valueLP);
			int curr = value - valueLP;
			if(curr < min) { min = curr; }
			if(curr > max) { max = curr; }
			maskImage.setPixelVal(i, j, curr);
		}
	}

	// normalize to range 0-255
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			maskImage.getPixelVal(i, j, value);
			int newvalue = 255l * (value - min) / (max - min);
			maskImage.setPixelVal(i, j, newvalue);
		}
	}
	
	// sharp mask file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_sharp_mask.pgm";
	char *newImageFile = new char[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, maskImage);
	delete[] newImageFile;
}
