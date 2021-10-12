#include "Correlation.h"

/**
 * This function computes for correlation of the image and the mask. It first takes in the different
 * mask sizes and image sizes to determine the offset and bounds of the window. It will make a window
 * and then multiply it by the mask's weights. This weighted sum is then returned.
 * @param: image ImageType reference, mask ImageType reference, integer values for original row and col
 * for current pixel, integer sizes for mask, integer sizes for original image
 * @return: integer weighted sum
 */
int computeCorrelation(ImageType& image, ImageType& mask, int row, int col, int Nm, int Mm, int N, int M) {
	// variables
	int k = row;
	int l = col;
	int om = Mm / 2;
	int on = Nm / 2;
	int value, valueM;
	int sum = 0;
	
	// grab sums for current window size and mask
	for(int i = 0; i < Mm; i++, k++) {
		l = col;
		for(int j = 0; j < Nm; j++, l++) {
			if((k-om>=0 && k-om<M) && (l-on>=0 && l-on<N)) {
				image.getPixelVal(k-om, l-on, value);
				mask.getPixelVal(i, j, valueM);
				sum += value * valueM;
			}
		}
	}
	
	return sum;
}

/**
 * This function takes in the image and mask's info and then calls computeCorrelation to compute
 * the correlation of the window and mask. It normalizes the values to 0-255 and then writes out
 * the image.
 * @param: fname character array to write image, image ImageType reference, mask ImageType
 * reference
 * @pre: all original values in variables
 * @post: correlated image written
 * @return: none
 */
void getCorrelation(char fname[], ImageType& image, ImageType& mask) {
	// variables
	int M, N, Q, Mm, Nm, value;
	image.getImageInfo(M, N, Q);
	mask.getImageInfo(Mm, Nm, Q);
	ImageType newImage(N, M, Q);

	// normalization map min max
	int min, max;
	image.getPixelVal(0, 0, min);
	max = min;
	
	// perform correlation
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < N; j++) {
			value = computeCorrelation(image, mask, i, j, Nm, Mm, N, M);
			if(value < min) { min = value; }
			if(value > max) { max = value; }
			newImage.setPixelVal(i, j, value);
		}
	}
	
	// normalize to range 0-255
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < N; j++) {
			newImage.getPixelVal(i, j, value);
			int newvalue = 255l * (value - min) / (max - min);
			newImage.setPixelVal(i, j, newvalue);
		}
	}

	// correlation file names and writing
	std::string newfname = "../images/" + std::string(fname) + "_correlation.pgm";
	char *newImageFile = new char[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
	delete[] newImageFile;
}
