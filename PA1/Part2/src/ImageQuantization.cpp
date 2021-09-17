#include "ImageQuantization.h"

/**
 * This function quantizes the number of gray levels in a source image and writes another
 * image for the quantization.
 * @param: fname character array to write image, image ImageType reference, bpp integer for
 * bits per pixel, mode integer to quantize at L-1 or 255.
 * @pre: image has data, quantized image written
 * @return: none
 */
void imageQuantization(char fname[], ImageType& image, int bpp, int mode) {
	// variables
	int M, N, Q, L, Q_new;
	int value;
	bool type;

	// original file names and reading
	std::string oldfname = "../images/" + std::string(fname) + ".pgm";
	char oldImageFile[oldfname.length() + 1];
	strcpy(oldImageFile, oldfname.c_str());
	readImageHeader(oldImageFile, N, M, Q, type);
	readImage(oldImageFile, image);

	// quantization properties
	L = pow(2, bpp);
	if(mode == 0) {
		// Q=L-1 (quantization reduced to bpp level)
		Q_new = L - 1;
	} else {
		// Q=255 (quantization keeping highest level)
		Q_new = 255;
	}
	ImageType newImage(N, M, Q_new);

	// perform quantization
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			int current = 0;
			image.getPixelVal(i, j, value);
			if(mode == 0) { current = value / pow(2, 8 - bpp); }
			else { current = value * pow(2, 8 - bpp); }
			newImage.setPixelVal(i, j, current);
		}
	}

	// quantized file names and writing
	std::string newfname = "../images/" + std::string(fname);
	if(mode == 0) { newfname += "_L" + std::to_string(L) + ".pgm"; }
	else { newfname += "_Q" + std::to_string(Q_new) + ".pgm"; }
	char newImageFile[newfname.length() + 1];
	strcpy(newImageFile, newfname.c_str());
	writeImage(newImageFile, newImage);
}
